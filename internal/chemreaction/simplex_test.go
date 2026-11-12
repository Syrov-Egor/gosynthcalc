package chemreaction

import (
	"context"
	"errors"
	"math/big"
	"slices"
	"testing"
)

func rationalMatrixFromInts(t *testing.T, rows [][]int64) [][]*big.Rat {
	t.Helper()
	if len(rows) == 0 {
		t.Fatal("rationalMatrixFromInts needs at least one row")
	}
	cols := len(rows[0])
	out := make([][]*big.Rat, len(rows))
	for i, row := range rows {
		if len(row) != cols {
			t.Fatalf("row %d has %d columns, want %d", i, len(row), cols)
		}
		out[i] = make([]*big.Rat, cols)
		for j, v := range row {
			out[i][j] = new(big.Rat).SetInt64(v)
		}
	}
	return out
}

func rationalVectorFromInts(vals ...int64) []*big.Rat {
	out := make([]*big.Rat, len(vals))
	for i, v := range vals {
		out[i] = new(big.Rat).SetInt64(v)
	}
	return out
}

func requireFeasibleSolution(t *testing.T, A [][]*big.Rat, b, c []*big.Rat) {
	t.Helper()
	if len(c) != len(A[0]) {
		t.Fatalf("solution length = %d, want %d", len(c), len(A[0]))
	}
	for j, v := range c {
		if v.Sign() < 0 {
			t.Fatalf("c[%d] = %s is negative", j, v)
		}
	}
	for i := range A {
		sum := new(big.Rat)
		for j := range c {
			sum.Add(sum, new(big.Rat).Mul(A[i][j], c[j]))
		}
		if sum.Cmp(b[i]) != 0 {
			t.Fatalf("row %d: A*c = %s, want %s", i, sum, b[i])
		}
	}
}

func TestPhaseIFeasible(t *testing.T) {
	A := rationalMatrixFromInts(t, [][]int64{{1, 1}})
	b := rationalVectorFromInts(2)
	c, feasible, err := phaseI(A, b)
	if err != nil {
		t.Fatal(err)
	}
	if !feasible {
		t.Fatal("phase I reported an infeasible system")
	}
	requireFeasibleSolution(t, A, b, c)
}

func TestPhaseIInfeasible(t *testing.T) {
	A := rationalMatrixFromInts(t, [][]int64{{1, 0}})
	b := rationalVectorFromInts(-1)
	_, feasible, err := phaseI(A, b)
	if err != nil {
		t.Fatal(err)
	}
	if feasible {
		t.Fatal("phase I reported a feasible system for x0 = -1 with x0 >= 0")
	}
}

func TestPhaseIZeroRow(t *testing.T) {
	A := rationalMatrixFromInts(t, [][]int64{{0, 0}})
	b := rationalVectorFromInts(0)
	c, feasible, err := phaseI(A, b)
	if err != nil {
		t.Fatal(err)
	}
	if !feasible {
		t.Fatal("a zero row with zero right-hand side must be feasible")
	}
	if !slices.EqualFunc(c, rationalVectorFromInts(0, 0), func(a, b *big.Rat) bool {
		return a.Cmp(b) == 0
	}) {
		t.Errorf("c = %v, want [0 0]", c)
	}
}

func TestPhaseIContextCancellation(t *testing.T) {
	ctx, cancel := context.WithCancel(context.Background())
	cancel()

	A := rationalMatrixFromInts(t, [][]int64{{1, 1}})
	b := rationalVectorFromInts(2)
	_, _, err := phaseI(A, b, simplexOptions{ctx: ctx})
	if !errors.Is(err, context.Canceled) {
		t.Errorf("phase I error = %v, want context.Canceled", err)
	}
}

func TestPhaseIIterationLimit(t *testing.T) {
	A := rationalMatrixFromInts(t, [][]int64{
		{1, 0, 1},
		{0, 1, 1},
	})
	b := rationalVectorFromInts(1, 1)
	_, _, err := phaseI(A, b, simplexOptions{maxIterations: 1})
	if !errors.Is(err, ErrSimplexIterationLimit) {
		t.Errorf("phase I error = %v, want ErrSimplexIterationLimit", err)
	}
}

func TestPhaseIMalformedInput(t *testing.T) {
	if _, _, err := phaseI(nil, nil); err != nil {
		t.Errorf("empty phase I error = %v, want nil", err)
	}
	A := rationalMatrixFromInts(t, [][]int64{{1, 2}})
	if _, _, err := phaseI(A, rationalVectorFromInts(1, 2)); !errors.Is(err, ErrInvalidMatrix) {
		t.Errorf("wrong RHS length error = %v, want ErrInvalidMatrix", err)
	}
	ragged := [][]*big.Rat{{big.NewRat(1, 1), big.NewRat(2, 1)}, {big.NewRat(3, 1)}}
	if _, _, err := phaseI(ragged, rationalVectorFromInts(1, 2)); !errors.Is(err, ErrInvalidMatrix) {
		t.Errorf("ragged matrix error = %v, want ErrInvalidMatrix", err)
	}
}

func TestBalanceExactSimpleReaction(t *testing.T) {
	reac, err := NewChemicalReaction("H2+O2=H2O")
	if err != nil {
		t.Fatal(err)
	}
	bal, err := reac.Balancer()
	if err != nil {
		t.Fatal(err)
	}
	A, _, err := bal.exactMatrixAndAnalysis()
	if err != nil {
		t.Fatal(err)
	}
	x, err := balanceExact(A)
	if err != nil {
		t.Fatal(err)
	}
	want := []*big.Rat{big.NewRat(2, 1), big.NewRat(1, 1), big.NewRat(2, 1)}
	if !slices.EqualFunc(x, want, func(a, b *big.Rat) bool { return a.Cmp(b) == 0 }) {
		t.Errorf("balanceExact = %v, want [2 1 2]", x)
	}
}

func TestBalanceExactEdgeCases(t *testing.T) {
	t.Run("no element rows", func(t *testing.T) {
		A := &exactMatrix{rows: 0, cols: 3, data: [][]*big.Rat{}}
		x, err := balanceExact(A)
		if err != nil {
			t.Fatal(err)
		}
		want := []*big.Rat{big.NewRat(1, 1), big.NewRat(1, 1), big.NewRat(1, 1)}
		if !slices.EqualFunc(x, want, func(a, b *big.Rat) bool { return a.Cmp(b) == 0 }) {
			t.Errorf("balanceExact = %v, want [1 1 1]", x)
		}
	})

	t.Run("zero row", func(t *testing.T) {
		A := exactMatrixFromInts(t, [][]int64{{0, 0, 0}})
		x, err := balanceExact(A)
		if err != nil {
			t.Fatal(err)
		}
		for i, v := range x {
			if v.Sign() <= 0 {
				t.Errorf("x[%d] = %s, want positive", i, v)
			}
		}
		if !exactBalancesRational(A, x) {
			t.Error("zero-row solution does not balance")
		}
	})

	t.Run("zero column is unconstrained", func(t *testing.T) {
		A := exactMatrixFromInts(t, [][]int64{{1, -1, 0}})
		x, err := balanceExact(A)
		if err != nil {
			t.Fatal(err)
		}
		for i, v := range x {
			if v.Sign() <= 0 {
				t.Errorf("x[%d] = %s, want positive", i, v)
			}
		}
		if !exactBalancesRational(A, x) {
			t.Error("zero-column solution does not balance")
		}
	})
}

func exactBalancesRational(A *exactMatrix, x []*big.Rat) bool {
	if A == nil || len(x) != A.cols {
		return false
	}
	for i := 0; i < A.rows; i++ {
		sum := new(big.Rat)
		for j := 0; j < A.cols; j++ {
			sum.Add(sum, new(big.Rat).Mul(A.data[i][j], x[j]))
		}
		if sum.Sign() != 0 {
			return false
		}
	}
	return true
}

func TestBalanceExactInfeasible(t *testing.T) {
	tests := []string{
		"H2+O2=H2O+Na",
		"FeS2+HNO3=Fe2(SO4)3+NO+H2SO4",
		"H2+O2=H2O+Na+Cl2+NaCl",
	}
	for _, reaction := range tests {
		t.Run(reaction, func(t *testing.T) {
			reac, err := NewChemicalReaction(reaction)
			if err != nil {
				t.Fatal(err)
			}
			bal, err := reac.Balancer()
			if err != nil {
				t.Fatal(err)
			}
			A, _, err := bal.exactMatrixAndAnalysis()
			if err != nil {
				t.Fatal(err)
			}
			if _, err := balanceExact(A); !errors.Is(err, ErrNoPositiveSolution) {
				t.Errorf("balanceExact error = %v, want ErrNoPositiveSolution", err)
			}
		})
	}
}

func TestBalanceExactContextCancellation(t *testing.T) {
	reac, err := NewChemicalReaction("KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O")
	if err != nil {
		t.Fatal(err)
	}
	bal, err := reac.Balancer()
	if err != nil {
		t.Fatal(err)
	}
	A, _, err := bal.exactMatrixAndAnalysis()
	if err != nil {
		t.Fatal(err)
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if _, err := balanceExact(A, simplexOptions{ctx: ctx}); !errors.Is(err, context.Canceled) {
		t.Errorf("balanceExact error = %v, want context.Canceled", err)
	}
}

func TestUniversalCorpus(t *testing.T) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		t.Fatal(err)
	}

	for i, rd := range reactions {
		reac, err := NewChemicalReaction(rd.reaction)
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		bal, err := reac.Balancer()
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		z, err := bal.UniversalExact(context.Background())
		if err != nil {
			t.Fatalf("row %d: UniversalExact: %v", i+1, err)
		}
		A, _, err := bal.exactMatrixAndAnalysis()
		if err != nil {
			t.Fatalf("row %d: analysis: %v", i+1, err)
		}
		requirePrimitivePositive(t, z, A)

		coefs, err := bal.Universal(context.Background())
		if err != nil {
			t.Fatalf("row %d: Universal: %v", i+1, err)
		}
		for j := range z {
			got := new(big.Rat).SetFloat64(coefs[j])
			want := new(big.Rat).SetInt(z[j])
			if got == nil || got.Cmp(want) != 0 {
				t.Errorf("row %d coefficient %d = %v, want %s", i+1, j, coefs[j], z[j])
			}
		}
	}
}

func TestUniversalExamples(t *testing.T) {
	tests := []struct {
		reaction string
		want     []int64
	}{
		{reaction: "H2+O2=H2O", want: []int64{2, 1, 2}},
		{reaction: "Bi2O3+TiO2=Bi4Ti3O12", want: []int64{2, 3, 1}},
		{reaction: "KMnO4+HCl=MnCl2+Cl2+H2O+KCl", want: []int64{2, 16, 2, 5, 8, 2}},
		{reaction: "H+O2=H1000001O2", want: []int64{1000001, 1, 1}},
	}
	for _, tt := range tests {
		t.Run(tt.reaction, func(t *testing.T) {
			reac, err := NewChemicalReaction(tt.reaction)
			if err != nil {
				t.Fatal(err)
			}
			bal, err := reac.Balancer()
			if err != nil {
				t.Fatal(err)
			}
			z, err := bal.UniversalExact(context.Background())
			if err != nil {
				t.Fatal(err)
			}
			want := make([]*big.Int, len(tt.want))
			for i, v := range tt.want {
				want[i] = big.NewInt(v)
			}
			if !slices.EqualFunc(z, want, func(a, b *big.Int) bool { return a.Cmp(b) == 0 }) {
				t.Errorf("UniversalExact = %v, want %v", z, tt.want)
			}
		})
	}
}

func TestUniversalHigherNullityExamples(t *testing.T) {
	tests := []struct {
		reaction string
		want     []int64
	}{
		{
			reaction: "KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O",
			want:     []int64{2, 2, 2, 1, 2, 1, 4},
		},
		{
			reaction: "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2",
			want:     []int64{4, 7, 1, 1, 1, 1, 5, 1},
		},
		{
			reaction: "Fe2O2H+FeO2H2=Fe3+Fe+Fe2O3H2",
			want:     []int64{8, 4, 1, 1, 8},
		},
	}
	for _, tt := range tests {
		t.Run(tt.reaction, func(t *testing.T) {
			reac, err := NewChemicalReaction(tt.reaction)
			if err != nil {
				t.Fatal(err)
			}
			bal, err := reac.Balancer()
			if err != nil {
				t.Fatal(err)
			}
			z, err := bal.UniversalExact(context.Background())
			if err != nil {
				t.Fatal(err)
			}
			A, a, err := bal.exactMatrixAndAnalysis()
			if err != nil {
				t.Fatal(err)
			}
			if a.nullity < 2 {
				t.Fatalf("nullity = %d, want at least 2", a.nullity)
			}
			requirePrimitivePositive(t, z, A)

			want := make([]*big.Int, len(tt.want))
			for i, v := range tt.want {
				want[i] = big.NewInt(v)
			}
			if !slices.EqualFunc(z, want, func(a, b *big.Int) bool { return a.Cmp(b) == 0 }) {
				t.Errorf("UniversalExact = %v, want %v", z, tt.want)
			}
		})
	}
}

func TestUniversalFractionalRows(t *testing.T) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		t.Fatal(err)
	}
	for i := 107; i <= 109; i++ {
		rd := reactions[i]
		t.Run(rd.reaction, func(t *testing.T) {
			reac, err := NewChemicalReaction(rd.reaction)
			if err != nil {
				t.Fatal(err)
			}
			bal, err := reac.Balancer()
			if err != nil {
				t.Fatal(err)
			}
			A, a, err := bal.exactMatrixAndAnalysis()
			if err != nil {
				t.Fatal(err)
			}
			if a.nullity != 1 {
				t.Fatalf("nullity = %d, want 1", a.nullity)
			}
			z, err := bal.UniversalExact(context.Background())
			if err != nil {
				t.Fatal(err)
			}
			requirePrimitivePositive(t, z, A)
			nullityZ, err := bal.solveNullityOneExact(a)
			if err != nil {
				t.Fatal(err)
			}
			if !slices.EqualFunc(z, nullityZ, func(a, b *big.Int) bool { return a.Cmp(b) == 0 }) {
				t.Errorf("exact LP = %v, nullity-one solver = %v", z, nullityZ)
			}
		})
	}
}

func TestUniversalContextCancellation(t *testing.T) {
	reac, err := NewChemicalReaction("KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O")
	if err != nil {
		t.Fatal(err)
	}
	bal, err := reac.Balancer()
	if err != nil {
		t.Fatal(err)
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if _, err := bal.UniversalExact(ctx); !errors.Is(err, context.Canceled) {
		t.Errorf("UniversalExact error = %v, want context.Canceled", err)
	}
}
