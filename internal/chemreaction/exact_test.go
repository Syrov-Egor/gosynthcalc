package chemreaction

import (
	"errors"
	"math"
	"math/big"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func exactMatrixFromInts(t *testing.T, rows [][]int64) *exactMatrix {
	t.Helper()
	if len(rows) == 0 {
		t.Fatal("exactMatrixFromInts needs at least one row")
	}
	cols := len(rows[0])
	data := make([][]*big.Rat, len(rows))
	for i, row := range rows {
		if len(row) != cols {
			t.Fatalf("row %d has %d columns, want %d", i, len(row), cols)
		}
		data[i] = make([]*big.Rat, cols)
		for j, v := range row {
			data[i][j] = new(big.Rat).SetInt64(v)
		}
	}
	return &exactMatrix{rows: len(rows), cols: cols, data: data}
}

func basisBalancesExactly(A *exactMatrix, basis [][]*big.Rat) bool {
	if len(basis) != A.cols {
		return false
	}
	if len(basis) == 0 {
		return true
	}
	k := len(basis[0])
	for j := 0; j < k; j++ {
		for i := 0; i < A.rows; i++ {
			sum := new(big.Rat)
			for c := 0; c < A.cols; c++ {
				sum.Add(sum, new(big.Rat).Mul(A.data[i][c], basis[c][j]))
			}
			if sum.Sign() != 0 {
				return false
			}
		}
	}
	return true
}

func TestRatFromFloat(t *testing.T) {
	tests := []struct {
		name string
		in   float64
		want string
	}{
		{name: "zero", in: 0, want: "0"},
		{name: "integer", in: 2, want: "2"},
		{name: "negative", in: -3.5, want: "-7/2"},
		{name: "decimal", in: 14.98, want: "749/50"},
		{name: "tenth", in: 0.1, want: "1/10"},
		{name: "negative zero", in: math.Copysign(0, -1), want: "0"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ratFromFloat(tt.in)
			if err != nil {
				t.Fatalf("ratFromFloat(%v) error: %v", tt.in, err)
			}
			want, _ := new(big.Rat).SetString(tt.want)
			if got.Cmp(want) != 0 {
				t.Errorf("ratFromFloat(%v) = %s, want %s", tt.in, got, want)
			}
		})
	}
}

func TestRatFromFloatRejectsNonFinite(t *testing.T) {
	for _, v := range []float64{math.NaN(), math.Inf(1), math.Inf(-1)} {
		if _, err := ratFromFloat(v); !errors.Is(err, ErrNonFiniteMatrix) {
			t.Errorf("ratFromFloat(%v) error = %v, want ErrNonFiniteMatrix", v, err)
		}
	}
}

func TestRatFromFloatDocumentsBinaryArtifact(t *testing.T) {
	// 0.1*3 is computed at run time as a float64, so the shortest decimal
	// representation is 0.30000000000000004, not 3/10. This is the known
	// nested-fractional-group limitation until the exact parser lands.
	x := 0.1
	got, err := ratFromFloat(x * 3)
	if err != nil {
		t.Fatal(err)
	}
	artifact, _ := new(big.Rat).SetString("0.30000000000000004")
	if got.Cmp(artifact) != 0 {
		t.Errorf("ratFromFloat(0.1*3) = %s, want the documented artifact %s", got, artifact)
	}
	if got.Cmp(big.NewRat(3, 10)) == 0 {
		t.Errorf("ratFromFloat(0.1*3) unexpectedly recovered 3/10")
	}
}

func TestRatFromDecimal(t *testing.T) {
	tests := []struct {
		in   string
		want string
	}{
		{in: "1.5", want: "3/2"},
		{in: "-0.25", want: "-1/4"},
		{in: "  12  ", want: "12"},
		{in: "7/3", want: "7/3"},
	}
	for _, tt := range tests {
		got, err := ratFromDecimal(tt.in)
		if err != nil {
			t.Errorf("ratFromDecimal(%q) error: %v", tt.in, err)
			continue
		}
		want, _ := new(big.Rat).SetString(tt.want)
		if got.Cmp(want) != 0 {
			t.Errorf("ratFromDecimal(%q) = %s, want %s", tt.in, got, want)
		}
	}

	if _, err := ratFromDecimal(""); !errors.Is(err, ErrInvalidMatrix) {
		t.Errorf("ratFromDecimal(\"\") error = %v, want ErrInvalidMatrix", err)
	}
	if _, err := ratFromDecimal("not a number"); !errors.Is(err, ErrInvalidMatrix) {
		t.Errorf("ratFromDecimal(bad) error = %v, want ErrInvalidMatrix", err)
	}
}

func TestBigIntsToFloat64(t *testing.T) {
	in := []*big.Int{
		big.NewInt(0),
		big.NewInt(1),
		big.NewInt(-5),
		new(big.Int).Lsh(big.NewInt(1), 40),
		new(big.Int).Lsh(big.NewInt(1), 53),
	}
	got, err := bigIntsToFloat64(in)
	if err != nil {
		t.Fatalf("bigIntsToFloat64 error: %v", err)
	}
	for i := range in {
		if new(big.Rat).SetInt(in[i]).Cmp(new(big.Rat).SetFloat64(got[i])) != 0 {
			t.Errorf("value %d was not converted exactly: %v", i, got[i])
		}
	}

	tooPrecise := new(big.Int).Add(new(big.Int).Lsh(big.NewInt(1), 53), big.NewInt(1))
	if _, err := bigIntsToFloat64([]*big.Int{tooPrecise}); !errors.Is(err, ErrCoefficientOverflow) {
		t.Errorf("2^53+1 error = %v, want ErrCoefficientOverflow", err)
	}
	huge := new(big.Int).Lsh(big.NewInt(1), 2000)
	if _, err := bigIntsToFloat64([]*big.Int{huge}); !errors.Is(err, ErrCoefficientOverflow) {
		t.Errorf("2^2000 error = %v, want ErrCoefficientOverflow", err)
	}
	if _, err := bigIntsToFloat64([]*big.Int{nil}); !errors.Is(err, ErrInvalidMatrix) {
		t.Errorf("nil coefficient error = %v, want ErrInvalidMatrix", err)
	}
}

func TestPrimitiveIntegers(t *testing.T) {
	got, err := primitiveIntegers([]*big.Rat{
		big.NewRat(1, 2),
		big.NewRat(1, 3),
		big.NewRat(1, 6),
	})
	if err != nil {
		t.Fatalf("primitiveIntegers error: %v", err)
	}
	want := []int64{3, 2, 1}
	for i := range want {
		if got[i].Cmp(big.NewInt(want[i])) != 0 {
			t.Errorf("primitiveIntegers[%d] = %s, want %d", i, got[i], want[i])
		}
	}

	mixed, err := primitiveIntegers([]*big.Rat{big.NewRat(1, 2), big.NewRat(2, 1)})
	if err != nil {
		t.Fatalf("primitiveIntegers error: %v", err)
	}
	if mixed[0].Cmp(big.NewInt(1)) != 0 || mixed[1].Cmp(big.NewInt(4)) != 0 {
		t.Errorf("primitiveIntegers([1/2,2]) = %v, want [1 4]", mixed)
	}

	if _, err := primitiveIntegers([]*big.Rat{big.NewRat(-1, 2)}); !errors.Is(err, ErrNoPositiveSolution) {
		t.Errorf("negative vector error = %v, want ErrNoPositiveSolution", err)
	}
	if _, err := primitiveIntegers([]*big.Rat{new(big.Rat)}); !errors.Is(err, ErrNoPositiveSolution) {
		t.Errorf("zero vector error = %v, want ErrNoPositiveSolution", err)
	}
	if _, err := primitiveIntegers(nil); !errors.Is(err, ErrInvalidMatrix) {
		t.Errorf("empty vector error = %v, want ErrInvalidMatrix", err)
	}
}

func TestExactBalances(t *testing.T) {
	A := exactMatrixFromInts(t, [][]int64{
		{2, 0, -2},
		{0, 2, -1},
	})
	if !exactBalances(A, []*big.Int{big.NewInt(2), big.NewInt(1), big.NewInt(2)}) {
		t.Error("2H2+O2=2H2O should balance exactly")
	}
	if exactBalances(A, []*big.Int{big.NewInt(1), big.NewInt(0), big.NewInt(1)}) {
		t.Error("H2+H2O=H2O is not balanced")
	}
	if exactBalances(A, []*big.Int{big.NewInt(2), big.NewInt(1)}) {
		t.Error("wrong-length vector must not balance")
	}
	if exactBalances(nil, []*big.Int{big.NewInt(1)}) {
		t.Error("nil matrix must not balance")
	}
}

func TestAnalyzeSmallMatrices(t *testing.T) {
	tests := []struct {
		name    string
		rows    [][]int64
		rank    int
		nullity int
	}{
		{
			name:    "zero matrix",
			rows:    [][]int64{{0, 0}, {0, 0}},
			rank:    0,
			nullity: 2,
		},
		{
			name:    "identity",
			rows:    [][]int64{{1, 0}, {0, 1}},
			rank:    2,
			nullity: 0,
		},
		{
			name:    "rank one line",
			rows:    [][]int64{{1, 2}, {2, 4}},
			rank:    1,
			nullity: 1,
		},
		{
			name: "thorne nullity two",
			rows: [][]int64{
				{1, 1, 0, -1, 0},
				{0, 0, 2, -4, -2},
				{1, 2, 0, 0, -1},
			},
			rank:    3,
			nullity: 2,
		},
		{
			name: "dependent element rows",
			rows: [][]int64{
				{1, 0, -1},
				{2, 0, -2},
			},
			rank:    1,
			nullity: 2,
		},
		{
			name: "duplicate columns",
			rows: [][]int64{
				{2, 2, -1, 0, 0},
				{6, 6, 0, -2, 0},
				{1, 1, 0, 0, -2},
			},
			rank:    3,
			nullity: 2,
		},
		{
			name: "square rank n-1",
			rows: [][]int64{
				{2, -2},
				{1, -1},
			},
			rank:    1,
			nullity: 1,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			A := exactMatrixFromInts(t, tt.rows)
			got, err := analyze(A)
			if err != nil {
				t.Fatalf("analyze error: %v", err)
			}
			if got.rank != tt.rank || got.nullity != tt.nullity {
				t.Errorf("rank/nullity = %d/%d, want %d/%d",
					got.rank, got.nullity, tt.rank, tt.nullity)
			}
			if got.rank+got.nullity != A.cols {
				t.Errorf("rank+nullity = %d, want %d", got.rank+got.nullity, A.cols)
			}
			if got.nullity == 0 {
				if got.basis != nil {
					t.Errorf("nullity 0 has a non-nil basis")
				}
				return
			}
			if !basisBalancesExactly(A, got.basis) {
				t.Errorf("null-space basis does not satisfy A*N = 0")
			}
		})
	}
}

func TestAnalyzeZeroMatrixBasisIsIdentity(t *testing.T) {
	A := exactMatrixFromInts(t, [][]int64{{0, 0, 0}})
	got, err := analyze(A)
	if err != nil {
		t.Fatal(err)
	}
	if got.nullity != 3 {
		t.Fatalf("nullity = %d, want 3", got.nullity)
	}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			want := int64(0)
			if i == j {
				want = 1
			}
			if got.basis[i][j].Cmp(new(big.Rat).SetInt64(want)) != 0 {
				t.Errorf("basis[%d][%d] = %s, want %d", i, j, got.basis[i][j], want)
			}
		}
	}
}

func TestNewExactMatrix(t *testing.T) {
	m := matDenseFromFloats(t, [][]float64{
		{2, 0, 2},
		{0, 2, 1},
	})
	A, err := newExactMatrix(m, 2)
	if err != nil {
		t.Fatalf("newExactMatrix error: %v", err)
	}
	if A.rows != 2 || A.cols != 3 {
		t.Fatalf("dims = %dx%d, want 2x3", A.rows, A.cols)
	}
	if A.data[0][0].Cmp(big.NewRat(2, 1)) != 0 {
		t.Errorf("reactant entry = %s, want 2", A.data[0][0])
	}
	if A.data[0][2].Cmp(big.NewRat(-2, 1)) != 0 {
		t.Errorf("product entry = %s, want -2", A.data[0][2])
	}

	nonFinite := matDenseFromFloats(t, [][]float64{{math.NaN()}})
	if _, err := newExactMatrix(nonFinite, 0); !errors.Is(err, ErrNonFiniteMatrix) {
		t.Errorf("NaN matrix error = %v, want ErrNonFiniteMatrix", err)
	}

	if _, err := newExactMatrix(nil, 0); !errors.Is(err, ErrInvalidMatrix) {
		t.Errorf("nil matrix error = %v, want ErrInvalidMatrix", err)
	}
	if _, err := newExactMatrix(m, 4); !errors.Is(err, ErrInvalidMatrix) {
		t.Errorf("bad separator error = %v, want ErrInvalidMatrix", err)
	}
}

func TestExactCorpusRankAndNullity(t *testing.T) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		t.Fatal(err)
	}
	if len(reactions) != 113 {
		t.Fatalf("corpus has %d reactions, want 113", len(reactions))
	}

	nullityCounts := map[int]int{}
	for i, rd := range reactions {
		reac, err := NewChemicalReaction(rd.reaction)
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		m, err := reac.Matrix()
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		A, err := newExactMatrix(m, reac.decomposer.separatorPos)
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		got, err := analyze(A)
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		if got.rank+got.nullity != A.cols {
			t.Errorf("row %d: rank+nullity = %d, want %d", i+1, got.rank+got.nullity, A.cols)
		}
		if got.nullity > 0 && !basisBalancesExactly(A, got.basis) {
			t.Errorf("row %d: basis does not satisfy A*N = 0", i+1)
		}
		if got.nullity == 1 {
			for _, r := range got.basisColumn(0) {
				if r.Sign() <= 0 {
					t.Errorf("row %d: nullity-1 direction is not strictly positive", i+1)
					break
				}
			}
		}
		nullityCounts[got.nullity]++
	}

	want := map[int]int{1: 103, 2: 7, 5: 3}
	for k, v := range want {
		if nullityCounts[k] != v {
			t.Errorf("nullity %d count = %d, want %d", k, nullityCounts[k], v)
		}
	}
	for k, v := range nullityCounts {
		if want[k] == 0 && v != 0 {
			t.Errorf("unexpected nullity %d count %d", k, v)
		}
	}
}

// basisColumn returns a copy of one column of the null-space basis.
func (a matrixAnalysis) basisColumn(j int) []*big.Rat {
	if j < 0 || len(a.basis) == 0 || j >= len(a.basis[0]) {
		return nil
	}
	out := make([]*big.Rat, len(a.basis))
	for i := range a.basis {
		out[i] = new(big.Rat).Set(a.basis[i][j])
	}
	return out
}

func matDenseFromFloats(t *testing.T, rows [][]float64) *mat.Dense {
	t.Helper()
	if len(rows) == 0 {
		t.Fatal("matDenseFromFloats needs at least one row")
	}
	cols := len(rows[0])
	for i := range rows {
		if len(rows[i]) != cols {
			t.Fatalf("row %d has %d columns, want %d", i, len(rows[i]), cols)
		}
	}
	m := mat.NewDense(len(rows), cols, nil)
	for i := range rows {
		for j := range rows[i] {
			m.Set(i, j, rows[i][j])
		}
	}
	return m
}
