package chemreaction

import (
	"errors"
	"math/big"
	"math/rand"
	"slices"
	"testing"
)

func TestPrimitiveScalingIdempotent(t *testing.T) {
	vectors := [][]*big.Rat{
		{big.NewRat(1, 2), big.NewRat(1, 3), big.NewRat(1, 6)},
		{big.NewRat(3, 1), big.NewRat(2, 1), big.NewRat(1, 1)},
		{big.NewRat(9999999999999999, 10000000000000000), big.NewRat(1, 1)},
	}
	for i, vector := range vectors {
		first, err := primitiveIntegers(vector)
		if err != nil {
			t.Fatalf("case %d: primitiveIntegers: %v", i, err)
		}
		rational := make([]*big.Rat, len(first))
		for j, v := range first {
			rational[j] = new(big.Rat).SetInt(v)
		}
		second, err := primitiveIntegers(rational)
		if err != nil {
			t.Fatalf("case %d: second primitiveIntegers: %v", i, err)
		}
		if !slices.EqualFunc(first, second, func(a, b *big.Int) bool { return a.Cmp(b) == 0 }) {
			t.Errorf("case %d: primitive scaling is not idempotent: %v vs %v", i, first, second)
		}
	}
}

func TestRandomFeasibleMatrices(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	for trial := 0; trial < 50; trial++ {
		rows := 1 + rng.Intn(4)
		cols := rows + 1 + rng.Intn(2)

		// Choose a known positive null vector. The last entry is fixed to 1
		// so the constructed matrix has integer entries.
		known := make([]*big.Rat, cols)
		coefficients := make([]int64, cols)
		for j := 0; j < cols-1; j++ {
			coefficients[j] = int64(1 + rng.Intn(5))
			known[j] = new(big.Rat).SetInt64(coefficients[j])
		}
		coefficients[cols-1] = 1
		known[cols-1] = big.NewRat(1, 1)

		data := make([][]int64, rows)
		for i := 0; i < rows; i++ {
			data[i] = make([]int64, cols)
			var weighted int64
			for j := 0; j < cols-1; j++ {
				data[i][j] = int64(rng.Intn(4))
				weighted += data[i][j] * coefficients[j]
			}
			data[i][cols-1] = -weighted
		}

		A := exactMatrixFromInts(t, data)
		if !exactBalancesRational(A, known) {
			t.Fatalf("trial %d: constructed matrix does not satisfy A*x = 0", trial)
		}

		x, err := balanceExact(A)
		if err != nil {
			t.Fatalf("trial %d: balanceExact: %v", trial, err)
		}
		if !exactBalancesRational(A, x) {
			t.Fatalf("trial %d: returned vector does not satisfy A*x = 0", trial)
		}
		for j, v := range x {
			if v.Sign() <= 0 {
				t.Fatalf("trial %d: x[%d] = %s, want positive", trial, j, v)
			}
		}
		z, err := primitiveIntegers(x)
		if err != nil {
			t.Fatalf("trial %d: primitiveIntegers: %v", trial, err)
		}
		if !exactBalances(A, z) {
			t.Fatalf("trial %d: primitive vector does not balance exactly", trial)
		}
	}
}

func TestRandomAllPositiveMatricesAreInfeasible(t *testing.T) {
	rng := rand.New(rand.NewSource(99))
	for trial := 0; trial < 30; trial++ {
		rows := 1 + rng.Intn(3)
		cols := rows + 1 + rng.Intn(3)
		data := make([][]int64, rows)
		for i := 0; i < rows; i++ {
			data[i] = make([]int64, cols)
			for j := 0; j < cols; j++ {
				data[i][j] = int64(1 + rng.Intn(5))
			}
		}
		A := exactMatrixFromInts(t, data)
		if _, err := balanceExact(A); !errors.Is(err, ErrNoPositiveSolution) {
			t.Fatalf("trial %d: all-positive matrix error = %v, want ErrNoPositiveSolution", trial, err)
		}
	}
}

func TestNormCoefficientsScaleInvariant(t *testing.T) {
	reac, err := NewChemicalReaction("H2+O2=H2O")
	if err != nil {
		t.Fatal(err)
	}

	if err := reac.SetCoefficients([]float64{2, 1, 2}); err != nil {
		t.Fatal(err)
	}
	first, err := reac.NormCoefficients()
	if err != nil {
		t.Fatal(err)
	}

	if err := reac.SetCoefficients([]float64{6, 3, 6}); err != nil {
		t.Fatal(err)
	}
	second, err := reac.NormCoefficients()
	if err != nil {
		t.Fatal(err)
	}

	if !slices.Equal(first, second) {
		t.Errorf("normalized coefficients changed under scaling: %v vs %v", first, second)
	}
}
