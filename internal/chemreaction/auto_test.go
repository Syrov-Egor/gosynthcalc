package chemreaction

import (
	"context"
	"errors"
	"math"
	"math/big"
	"slices"
	"testing"
)

func exactIntsFromFloatCoefs(t *testing.T, coefs []float64) []*big.Int {
	t.Helper()
	out := make([]*big.Int, len(coefs))
	for i, f := range coefs {
		if math.IsNaN(f) || math.IsInf(f, 0) {
			t.Fatalf("coefficient %d is not finite: %v", i, f)
		}
		r := new(big.Rat).SetFloat64(f)
		if r == nil || !r.IsInt() {
			t.Fatalf("coefficient %d = %v is not an integer", i, f)
		}
		out[i] = new(big.Int).Set(r.Num())
	}
	return out
}

func TestAutoCorpusExactRouting(t *testing.T) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		t.Fatal(err)
	}

	methodCounts := map[string]int{}
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

		result, err := bal.Auto()
		if err != nil {
			t.Fatalf("row %d: Auto: %v", i+1, err)
		}
		wantMethod := MethodExactLP
		if analysis.nullity == 1 {
			wantMethod = MethodInverse
		}
		if result.Method != wantMethod {
			t.Errorf("row %d: method = %q, want %q", i+1, result.Method, wantMethod)
		}
		methodCounts[result.Method]++

		z := exactIntsFromFloatCoefs(t, result.Result)
		requirePrimitivePositive(t, z, A)
	}

	if methodCounts[MethodInverse] != 103 {
		t.Errorf("inverse method count = %d, want 103", methodCounts[MethodInverse])
	}
	if methodCounts[MethodExactLP] != 10 {
		t.Errorf("exact LP method count = %d, want 10", methodCounts[MethodExactLP])
	}
}

func TestAutoDeterministic(t *testing.T) {
	const reaction = "KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O"
	first, err := NewChemicalReaction(reaction)
	if err != nil {
		t.Fatal(err)
	}
	firstBal, err := first.Balancer()
	if err != nil {
		t.Fatal(err)
	}
	second, err := NewChemicalReaction(reaction)
	if err != nil {
		t.Fatal(err)
	}
	secondBal, err := second.Balancer()
	if err != nil {
		t.Fatal(err)
	}

	firstResult, err := firstBal.Auto()
	if err != nil {
		t.Fatal(err)
	}
	secondResult, err := secondBal.Auto()
	if err != nil {
		t.Fatal(err)
	}
	if firstResult.Method != secondResult.Method {
		t.Errorf("method differs: %q vs %q", firstResult.Method, secondResult.Method)
	}
	if !slices.Equal(firstResult.Result, secondResult.Result) {
		t.Errorf("result differs: %v vs %v", firstResult.Result, secondResult.Result)
	}
}

func TestAutoInfeasible(t *testing.T) {
	reactions := []string{
		"H2+O2=H2O+Na",
		"H2+H2O=H2O2",
		"FeS2+HNO3=Fe2(SO4)3+NO+H2SO4",
		"H2+O2=H2O+Na+Cl2+NaCl",
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
			_, err = bal.Auto()
			if !errors.Is(err, ErrNoPositiveSolution) {
				t.Errorf("Auto error = %v, want ErrNoPositiveSolution", err)
			}
		})
	}
}

func TestAutoNullityZero(t *testing.T) {
	m := matDenseFromFloats(t, [][]float64{
		{2, 2},
		{1, 2},
	})
	bal := newBalancer(m, 1, true, 8, StrategyExact, 1e-8)
	_, err := bal.Auto()
	if !errors.Is(err, ErrNoPositiveSolution) {
		t.Errorf("Auto error = %v, want ErrNoPositiveSolution", err)
	}
}

func TestAutoContextCancellation(t *testing.T) {
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
	if _, err := bal.AutoContext(ctx); !errors.Is(err, context.Canceled) {
		t.Errorf("AutoContext error = %v, want context.Canceled", err)
	}
}

func TestStrategyLegacyCorpus(t *testing.T) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		t.Fatal(err)
	}

	options := ReacOptions{
		Rmode:      Balance,
		Target:     DefaultTarget,
		TargetMass: DefaultTargetMass,
		Intify:     DefaultIntify,
		Precision:  DefaultPrecision,
		Tolerance:  DefaultTolerance,
		Strategy:   StrategyLegacy,
	}
	for i, rd := range reactions {
		reac, err := NewChemicalReaction(rd.reaction, options)
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		bal, err := reac.Balancer()
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		result, err := bal.Auto()
		if err != nil {
			t.Fatalf("row %d: legacy Auto: %v", i+1, err)
		}
		if !bal.validFloatResult(result.Result) {
			t.Errorf("row %d: legacy result %v does not validate", i+1, result.Result)
		}
	}
}

func TestStrategyLegacyFallsBackToExact(t *testing.T) {
	// All three historical heuristics fail on this nullity-2 reaction, so the
	// legacy strategy must still succeed through the exact simplex.
	const reaction = "Fe2O2H+FeO2H2=Fe3+Fe+Fe2O3H2"
	reac, err := NewChemicalReaction(reaction)
	if err != nil {
		t.Fatal(err)
	}
	bal, err := reac.Balancer()
	if err != nil {
		t.Fatal(err)
	}
	bal.strategy = StrategyLegacy

	result, err := bal.Auto()
	if err != nil {
		t.Fatal(err)
	}
	if result.Method != MethodExactLP {
		t.Errorf("method = %q, want %q", result.Method, MethodExactLP)
	}
	A, _, err := bal.exactMatrixAndAnalysis()
	if err != nil {
		t.Fatal(err)
	}
	z := exactIntsFromFloatCoefs(t, result.Result)
	requirePrimitivePositive(t, z, A)
}

func TestAutoExactIgnoresLegacyControls(t *testing.T) {
	const reaction = "KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O"

	strict := ReacOptions{
		Rmode:      Balance,
		Target:     DefaultTarget,
		TargetMass: DefaultTargetMass,
		Intify:     true,
		Precision:  2,
		Tolerance:  1e-12,
		Strategy:   StrategyExact,
	}
	loose := ReacOptions{
		Rmode:      Balance,
		Target:     DefaultTarget,
		TargetMass: DefaultTargetMass,
		Intify:     false,
		Precision:  8,
		Tolerance:  1e-3,
		Strategy:   StrategyExact,
	}

	first, err := NewChemicalReaction(reaction, strict)
	if err != nil {
		t.Fatal(err)
	}
	second, err := NewChemicalReaction(reaction, loose)
	if err != nil {
		t.Fatal(err)
	}
	firstResult, err := first.Coefficients()
	if err != nil {
		t.Fatal(err)
	}
	secondResult, err := second.Coefficients()
	if err != nil {
		t.Fatal(err)
	}

	if firstResult.Method != MethodExactLP || secondResult.Method != MethodExactLP {
		t.Fatalf("methods = %q/%q, want %q",
			firstResult.Method, secondResult.Method, MethodExactLP)
	}
	if !slices.Equal(firstResult.Result, secondResult.Result) {
		t.Errorf("exact results depend on Intify/Precision/Tolerance: %v vs %v",
			firstResult.Result, secondResult.Result)
	}
	for i, coef := range firstResult.Result {
		if coef != math.Trunc(coef) {
			t.Errorf("coefficient %d = %v is not an integer", i, coef)
		}
	}
}
