package chemreaction

import (
	"errors"
	"math/big"
	"slices"
	"testing"
)

func bigIntsEqual(got []*big.Int, want []int64) bool {
	if len(got) != len(want) {
		return false
	}
	for i := range got {
		if got[i].Cmp(big.NewInt(want[i])) != 0 {
			return false
		}
	}
	return true
}

func TestIsBalancedExactIntegerValidation(t *testing.T) {
	// Check mode with a very large tolerance would accept [1,1,1] as
	// balanced by float comparison. IsBalanced sees integer coefficients and
	// an exact matrix, so it must reject them exactly.
	options := ReacOptions{
		Rmode:      Check,
		Target:     DefaultTarget,
		TargetMass: DefaultTargetMass,
		Intify:     true,
		Precision:  DefaultPrecision,
		Tolerance:  10,
		Strategy:   StrategyExact,
	}
	reac, err := NewChemicalReaction("H2+O2=H2O", options)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := reac.Coefficients(); err != nil {
		t.Fatalf("Check mode should accept the coefficients under the large scaled tolerance: %v", err)
	}
	if reac.IsBalanced() {
		t.Error("IsBalanced() = true for integer coefficients that do not balance exactly")
	}
}

func TestIsBalancedFractionalFallback(t *testing.T) {
	options := ReacOptions{
		Rmode:      Check,
		Target:     DefaultTarget,
		TargetMass: DefaultTargetMass,
		Intify:     true,
		Precision:  DefaultPrecision,
		Tolerance:  DefaultTolerance,
		Strategy:   StrategyExact,
	}
	reac, err := NewChemicalReaction("0.5H2+0.25O2=0.5H2O", options)
	if err != nil {
		t.Fatal(err)
	}
	if !reac.IsBalanced() {
		t.Error("IsBalanced() = false for exactly balanced fractional coefficients")
	}
}

func TestCoefficientsExact(t *testing.T) {
	tests := []struct {
		reaction string
		method   string
		want     []int64
	}{
		{
			reaction: "H2+O2=H2O",
			method:   MethodInverse,
			want:     []int64{2, 1, 2},
		},
		{
			reaction: "KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O",
			method:   MethodExactLP,
			want:     []int64{2, 2, 2, 1, 2, 1, 4},
		},
		{
			reaction: "H+O2=H1000001O2",
			method:   MethodInverse,
			want:     []int64{1000001, 1, 1},
		},
	}
	for _, tt := range tests {
		t.Run(tt.reaction, func(t *testing.T) {
			reac, err := NewChemicalReaction(tt.reaction)
			if err != nil {
				t.Fatal(err)
			}
			got, err := reac.CoefficientsExact()
			if err != nil {
				t.Fatalf("CoefficientsExact() error: %v", err)
			}
			if got.Method != tt.method {
				t.Errorf("method = %q, want %q", got.Method, tt.method)
			}
			if !bigIntsEqual(got.Result, tt.want) {
				t.Errorf("result = %v, want %v", got.Result, tt.want)
			}
		})
	}
}

func TestCoefficientsExactAfterFloatAccess(t *testing.T) {
	reac, err := NewChemicalReaction("KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O")
	if err != nil {
		t.Fatal(err)
	}
	if _, err := reac.Coefficients(); err != nil {
		t.Fatal(err)
	}
	got, err := reac.CoefficientsExact()
	if err != nil {
		t.Fatal(err)
	}
	want := []int64{2, 2, 2, 1, 2, 1, 4}
	if !bigIntsEqual(got.Result, want) {
		t.Errorf("result = %v, want %v", got.Result, want)
	}
	if got.Method != MethodExactLP {
		t.Errorf("method = %q, want %q", got.Method, MethodExactLP)
	}
}

func TestCoefficientsExactUserCoefficients(t *testing.T) {
	options := ReacOptions{
		Rmode:      Check,
		Target:     DefaultTarget,
		TargetMass: DefaultTargetMass,
		Intify:     true,
		Precision:  DefaultPrecision,
		Tolerance:  DefaultTolerance,
		Strategy:   StrategyExact,
	}
	reac, err := NewChemicalReaction("0.5H2+0.25O2=0.5H2O", options)
	if err != nil {
		t.Fatal(err)
	}
	got, err := reac.CoefficientsExact()
	if err != nil {
		t.Fatal(err)
	}
	if got.Method != MethodUser {
		t.Errorf("method = %q, want %q", got.Method, MethodUser)
	}
	if !bigIntsEqual(got.Result, []int64{2, 1, 2}) {
		t.Errorf("result = %v, want [2 1 2]", got.Result)
	}
}

func TestCoefficientsExactAvoidsFloatOverflow(t *testing.T) {
	// The exact solution contains 9999999999999999, an odd integer above
	// 2^53 that cannot be represented exactly as a float64. The float API
	// must report overflow instead of a corrupted coefficient, while the
	// exact accessor must return the primitive integers.
	reac, err := NewChemicalReaction("H+O2=H0.9999999999999999O2")
	if err != nil {
		t.Fatal(err)
	}

	if _, err := reac.Coefficients(); !errors.Is(err, ErrCoefficientOverflow) {
		t.Errorf("Coefficients() error = %v, want ErrCoefficientOverflow", err)
	}

	got, err := reac.CoefficientsExact()
	if err != nil {
		t.Fatalf("CoefficientsExact() error: %v", err)
	}
	want := []int64{9999999999999999, 10000000000000000, 10000000000000000}
	if !bigIntsEqual(got.Result, want) {
		t.Errorf("CoefficientsExact() = %v, want %v", got.Result, want)
	}
	if got.Method != MethodInverse {
		t.Errorf("method = %q, want %q", got.Method, MethodInverse)
	}
}

func TestCoefficientsExactInfeasible(t *testing.T) {
	reac, err := NewChemicalReaction("FeS2+HNO3=Fe2(SO4)3+NO+H2SO4")
	if err != nil {
		t.Fatal(err)
	}
	_, err = reac.CoefficientsExact()
	if !errors.Is(err, ErrNoPositiveSolution) {
		t.Errorf("CoefficientsExact() error = %v, want ErrNoPositiveSolution", err)
	}
}

func TestErrorClassificationIsDistinct(t *testing.T) {
	if errors.Is(ErrMethodNotPositive, ErrNoPositiveSolution) {
		t.Error("a failed heuristic must not be reported as infeasibility")
	}
	if errors.Is(ErrCoefficientOverflow, ErrNoPositiveSolution) {
		t.Error("overflow must not be reported as infeasibility")
	}
}

func TestNormCoefficientsFromExactResult(t *testing.T) {
	reac, err := NewChemicalReaction("KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O")
	if err != nil {
		t.Fatal(err)
	}
	norm, err := reac.NormCoefficients()
	if err != nil {
		t.Fatal(err)
	}
	if len(norm) != 7 {
		t.Fatalf("normalized coefficients length = %d, want 7", len(norm))
	}
	if norm[0] != 2 {
		t.Errorf("normalized first coefficient = %v, want 2", norm[0])
	}
	if !slices.Equal(norm, []float64{2, 2, 2, 1, 2, 1, 4}) {
		t.Errorf("NormCoefficients() = %v, want [2 2 2 1 2 1 4]", norm)
	}
}
