package chemreaction

import (
	"errors"
	"math"
	"math/big"
	"slices"
	"testing"
)

func requirePrimitivePositive(t *testing.T, z []*big.Int, A *exactMatrix) {
	t.Helper()
	if len(z) != A.cols {
		t.Fatalf("coefficient count = %d, want %d", len(z), A.cols)
	}
	gcd := new(big.Int)
	for i, v := range z {
		if v == nil || v.Sign() <= 0 {
			t.Fatalf("coefficient %d = %v, want a positive integer", i, v)
		}
		gcd.GCD(nil, nil, gcd, new(big.Int).Abs(v))
	}
	if gcd.Cmp(big.NewInt(1)) != 0 {
		t.Errorf("coefficients are not primitive: gcd = %s", gcd)
	}
	if !exactBalances(A, z) {
		t.Errorf("coefficients %v do not satisfy A*z = 0 exactly", z)
	}
}

func TestNullityOneSolution(t *testing.T) {
	rat := func(n, d int64) *big.Rat { return big.NewRat(n, d) }

	t.Run("positive direction", func(t *testing.T) {
		a := matrixAnalysis{
			nullity: 1,
			basis:   [][]*big.Rat{{rat(1, 2)}, {rat(1, 3)}},
		}
		got, err := nullityOneSolution(a)
		if err != nil {
			t.Fatal(err)
		}
		if got[0].Cmp(rat(1, 2)) != 0 || got[1].Cmp(rat(1, 3)) != 0 {
			t.Errorf("got %v, want the original positive direction", got)
		}
	})

	t.Run("negative direction is flipped", func(t *testing.T) {
		a := matrixAnalysis{
			nullity: 1,
			basis:   [][]*big.Rat{{rat(-1, 2)}, {rat(-1, 3)}},
		}
		got, err := nullityOneSolution(a)
		if err != nil {
			t.Fatal(err)
		}
		if got[0].Cmp(rat(1, 2)) != 0 || got[1].Cmp(rat(1, 3)) != 0 {
			t.Errorf("got %v, want the negated direction", got)
		}
	})

	t.Run("zero entry", func(t *testing.T) {
		a := matrixAnalysis{
			nullity: 1,
			basis:   [][]*big.Rat{{rat(2, 1)}, {new(big.Rat)}},
		}
		if _, err := nullityOneSolution(a); !errors.Is(err, ErrNoPositiveSolution) {
			t.Errorf("error = %v, want ErrNoPositiveSolution", err)
		}
	})

	t.Run("mixed signs", func(t *testing.T) {
		a := matrixAnalysis{
			nullity: 1,
			basis:   [][]*big.Rat{{rat(2, 1)}, {rat(-1, 1)}},
		}
		if _, err := nullityOneSolution(a); !errors.Is(err, ErrNoPositiveSolution) {
			t.Errorf("error = %v, want ErrNoPositiveSolution", err)
		}
	})

	t.Run("nullity zero", func(t *testing.T) {
		if _, err := nullityOneSolution(matrixAnalysis{nullity: 0}); !errors.Is(err, ErrNoPositiveSolution) {
			t.Errorf("error = %v, want ErrNoPositiveSolution", err)
		}
	})

	t.Run("nullity two", func(t *testing.T) {
		if _, err := nullityOneSolution(matrixAnalysis{nullity: 2}); !errors.Is(err, ErrNullityNotOne) {
			t.Errorf("error = %v, want ErrNullityNotOne", err)
		}
	})

	t.Run("malformed basis", func(t *testing.T) {
		if _, err := nullityOneSolution(matrixAnalysis{nullity: 1}); !errors.Is(err, ErrInvalidMatrix) {
			t.Errorf("error = %v, want ErrInvalidMatrix", err)
		}
	})
}

func TestInvSimpleReactions(t *testing.T) {
	tests := []struct {
		reaction string
		want     []float64
	}{
		{reaction: "H2+O2=H2O", want: []float64{2, 1, 2}},
		{reaction: "Bi2O3+TiO2=Bi4Ti3O12", want: []float64{2, 3, 1}},
		{reaction: "KMnO4+HCl=MnCl2+Cl2+H2O+KCl", want: []float64{2, 16, 2, 5, 8, 2}},
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
			got, err := bal.Inv()
			if err != nil {
				t.Fatalf("Inv() error: %v", err)
			}
			if !slices.Equal(got, tt.want) {
				t.Errorf("Inv() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestInvRejectsNonPositiveNullityOne(t *testing.T) {
	tests := []struct {
		name     string
		reaction string
	}{
		{name: "zero coefficient", reaction: "H2+O2=H2O+Na"},
		{name: "mixed signs", reaction: "H2+H2O=H2O2"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			reac, err := NewChemicalReaction(tt.reaction)
			if err != nil {
				t.Fatal(err)
			}
			bal, err := reac.Balancer()
			if err != nil {
				t.Fatal(err)
			}
			_, err = bal.Inv()
			if !errors.Is(err, ErrNoPositiveSolution) {
				t.Errorf("Inv() error = %v, want ErrNoPositiveSolution", err)
			}
		})
	}
}

func TestInvSquareRankNMinusOne(t *testing.T) {
	// The unsigned matrix of H2O = H2O: two elements, two identical
	// composition columns. Signed and with separator 1 it has rank 1 and
	// nullity 1, and the unique null direction is [1, 1].
	m := matDenseFromFloats(t, [][]float64{
		{2, 2},
		{1, 1},
	})
	bal := newBalancer(m, 1, true, 8, StrategyExact, 1e-8)
	got, err := bal.Inv()
	if err != nil {
		t.Fatalf("Inv() error: %v", err)
	}
	want := []float64{1, 1}
	if !slices.Equal(got, want) {
		t.Errorf("Inv() = %v, want %v", got, want)
	}
}

func TestInvRejectsNullityZero(t *testing.T) {
	// The unsigned matrix of H2O = H2O2: signed with separator 1 it is
	// [[2,-2],[1,-2]], which has full rank, so there is no non-zero null
	// vector at all.
	m := matDenseFromFloats(t, [][]float64{
		{2, 2},
		{1, 2},
	})
	bal := newBalancer(m, 1, true, 8, StrategyExact, 1e-8)
	_, err := bal.Inv()
	if !errors.Is(err, ErrNoPositiveSolution) {
		t.Errorf("Inv() error = %v, want ErrNoPositiveSolution", err)
	}
}

func TestInvRejectsHigherNullity(t *testing.T) {
	reactions := []string{
		"KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O",
		"C2H6O+C2H6O=C+H2+O2",
		"H2+O2=H2O+Na+Cl2+NaCl",
		"Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2",
	}
	for _, reaction := range reactions {
		t.Run(reaction, func(t *testing.T) {
			reac, err := NewChemicalReaction(reaction)
			if err != nil {
				t.Fatal(err)
			}
			bal, err := reac.Balancer()
			if err != nil {
				t.Fatal(err)
			}
			_, err = bal.Inv()
			if !errors.Is(err, ErrNullityNotOne) {
				t.Errorf("Inv() error = %v, want ErrNullityNotOne", err)
			}
		})
	}
}

func TestInvAllNullityOneCorpusRows(t *testing.T) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		t.Fatal(err)
	}

	checked := 0
	for i, rd := range reactions {
		reac, err := NewChemicalReaction(rd.reaction)
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		bal, err := reac.Balancer()
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		A, analysis, err := bal.exactMatrixAndAnalysis()
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		if analysis.nullity != 1 {
			continue
		}
		checked++

		z, err := bal.solveNullityOneExact(analysis)
		if err != nil {
			t.Fatalf("row %d: solveNullityOneExact: %v", i+1, err)
		}
		requirePrimitivePositive(t, z, A)

		coefs, err := bal.Inv()
		if err != nil {
			t.Fatalf("row %d: Inv: %v", i+1, err)
		}
		if len(coefs) != len(z) {
			t.Fatalf("row %d: Inv returned %d coefficients, want %d", i+1, len(coefs), len(z))
		}
		for j := range z {
			got := new(big.Rat).SetFloat64(coefs[j])
			want := new(big.Rat).SetInt(z[j])
			if got == nil || got.Cmp(want) != 0 {
				t.Errorf("row %d coefficient %d = %v is not the exact primitive value %s",
					i+1, j, coefs[j], z[j])
			}
		}
	}
	if checked != 103 {
		t.Errorf("checked %d nullity-1 rows, want 103", checked)
	}
}

func TestInvFractionalRowsUseExactIntegers(t *testing.T) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		t.Fatal(err)
	}

	// Rows 108-110 (0-based 107-109) carry decimal atom counts such as
	// O14.98, Br1.99 and K1.97.
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
			A, analysis, err := bal.exactMatrixAndAnalysis()
			if err != nil {
				t.Fatal(err)
			}
			if analysis.nullity != 1 {
				t.Fatalf("nullity = %d, want 1", analysis.nullity)
			}
			z, err := bal.solveNullityOneExact(analysis)
			if err != nil {
				t.Fatalf("solveNullityOneExact error: %v", err)
			}
			requirePrimitivePositive(t, z, A)

			coefs, err := bal.Inv()
			if err != nil {
				t.Fatalf("Inv error: %v", err)
			}
			if slices.Equal(coefs, rd.coefs) {
				t.Errorf("Inv returned the old heuristic floats %v, want exact primitive integers", coefs)
			}
			for j, coef := range coefs {
				if coef != math.Trunc(coef) {
					t.Errorf("coefficient %d = %v is not an integer", j, coef)
				}
			}
		})
	}
}
