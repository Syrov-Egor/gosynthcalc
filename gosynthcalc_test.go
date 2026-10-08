package gosynthcalc

import (
	"errors"
	"fmt"
	"io"
	"os"
	"slices"
	"strings"
	"testing"
)

type benchmarkSkippedReaction struct {
	line     int
	reaction string
	err      error
}

type benchmarkData struct {
	formulas  []string
	reactions []string
	skipped   []benchmarkSkippedReaction
}

func loadBenchmarkData(reader io.Reader) (benchmarkData, error) {
	data := benchmarkData{}
	b, err := io.ReadAll(reader)
	if err != nil {
		return data, fmt.Errorf("read benchmark dataset: %w", err)
	}
	reactionsStr := strings.Split(string(b), "\n")

	for i, reac := range reactionsStr {
		if strings.TrimSpace(reac) == "" {
			continue
		}
		reacO, err := NewChemicalReaction(reac)
		if err != nil {
			data.skipped = append(data.skipped, benchmarkSkippedReaction{i + 1, reac, err})
			continue
		}

		forms, err := reacO.ChemFormulas()
		if err != nil {
			data.skipped = append(data.skipped, benchmarkSkippedReaction{i + 1, reac, err})
			continue
		}
		data.reactions = append(data.reactions, reac)
		for _, f := range forms {
			data.formulas = append(data.formulas, f.Formula())
		}
	}
	if len(data.reactions) == 0 || len(data.formulas) == 0 {
		return data, fmt.Errorf("no valid benchmark reactions")
	}
	return data, nil
}

func setup(tb testing.TB, fname string) ([]string, []string) {
	tb.Helper()
	file, err := os.Open(fname)
	if err != nil {
		tb.Fatalf("open benchmark dataset %q: %v", fname, err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			tb.Errorf("close benchmark dataset %q: %v", fname, err)
		}
	}()

	data, err := loadBenchmarkData(file)
	for _, skipped := range data.skipped[:min(len(data.skipped), 5)] {
		tb.Logf("%s:%d: skipped reaction %q: %v", fname, skipped.line, skipped.reaction, skipped.err)
	}
	if err != nil {
		tb.Fatalf("load benchmark dataset %q: %v", fname, err)
	}
	return data.formulas, data.reactions
}

func BenchmarkChemicalFormula_output(b *testing.B) {
	b.ReportAllocs()
	formulas, _ := setup(b, "data/text_mined_reactions.txt")
	n := len(formulas)

	i := 0
	for b.Loop() {
		form := formulas[i%n]
		i++

		formulaObj, err := NewChemicalFormula(form)
		if err != nil {
			b.Fatalf("formula %q: %v", form, err)
		}
		_ = formulaObj.Output()
	}
}

func BenchmarkChemicalReaction_output(b *testing.B) {
	b.ReportAllocs()
	_, reactions := setup(b, "data/text_mined_reactions.txt")
	n := len(reactions)

	i := 0
	for b.Loop() {
		reac := reactions[i%n]
		i++

		reactionObj, err := NewChemicalReaction(reac)
		if err != nil {
			b.Fatalf("reaction %q: %v", reac, err)
		}
		if _, err := reactionObj.Output(); err != nil {
			b.Fatalf("reaction %q: %v", reac, err)
		}
	}
}

func TestLoadBenchmarkData(t *testing.T) {
	input := "\n2H2+O2==2H2O\nH2+O2==H1.2.3O\nBaCuO2(011)+Y2O3==Y2BaCuO5\nnot a reaction\n\nNaOH+HCl==NaCl+H2O\n"
	data, err := loadBenchmarkData(strings.NewReader(input))
	if err != nil {
		t.Fatalf("loadBenchmarkData() unexpected error: %v", err)
	}
	wantReactions := []string{"2H2+O2==2H2O", "NaOH+HCl==NaCl+H2O"}
	wantFormulas := []string{"H2", "O2", "H2O", "NaOH", "HCl", "NaCl", "H2O"}
	if !slices.Equal(data.reactions, wantReactions) {
		t.Errorf("reactions = %v, expected %v", data.reactions, wantReactions)
	}
	if !slices.Equal(data.formulas, wantFormulas) {
		t.Errorf("formulas = %v, expected %v", data.formulas, wantFormulas)
	}
	if len(data.skipped) != 3 {
		t.Fatalf("skipped %d records, expected 3", len(data.skipped))
	}
	for i, skipped := range data.skipped {
		if skipped.line != i+3 || skipped.reaction == "" || skipped.err == nil {
			t.Errorf("incomplete skipped-record diagnostic: %+v", skipped)
		}
	}
}

func TestLoadBenchmarkData_empty(t *testing.T) {
	for _, input := range []string{"", " \n\n", "H2+O2==H1.2.3O\nnot a reaction\n"} {
		t.Run(input, func(t *testing.T) {
			data, err := loadBenchmarkData(strings.NewReader(input))
			if err == nil || !strings.Contains(err.Error(), "no valid benchmark reactions") {
				t.Errorf("loadBenchmarkData() error = %v, expected no valid benchmark reactions", err)
			}
			if len(data.formulas) != 0 || len(data.reactions) != 0 {
				t.Errorf("loadBenchmarkData() retained invalid records: %+v", data)
			}
		})
	}
}

type benchmarkErrorReader struct{}

func (benchmarkErrorReader) Read([]byte) (int, error) {
	return 0, io.ErrUnexpectedEOF
}

func TestLoadBenchmarkData_readError(t *testing.T) {
	_, err := loadBenchmarkData(benchmarkErrorReader{})
	if !errors.Is(err, io.ErrUnexpectedEOF) {
		t.Errorf("loadBenchmarkData() error = %v, expected unexpected EOF", err)
	}
}

func TestBenchmarkDataset(t *testing.T) {
	// Exercise the same corpus-loading path as benchmarks in normal test runs.
	formulas, reactions := setup(t, "data/text_mined_reactions.txt")
	if len(formulas) == 0 || len(reactions) == 0 {
		t.Fatal("benchmark dataset has no usable records")
	}
}

func ExampleNewChemicalFormula() {
	formula, err := NewChemicalFormula("H2O")
	if err != nil {
		fmt.Printf("Error: %v\n", err)
		return
	}
	molarMass := formula.MolarMass()
	fmt.Printf("%v", molarMass)
	// Output: 18.015
}

func ExampleNewChemicalReaction() {
	reaction, err := NewChemicalReaction("FeCl3+SO2+H2O=FeCl2+H2SO4+HCl")
	if err != nil {
		fmt.Printf("Error: %v\n", err)
		return
	}
	coefs, err := reaction.Coefficients()
	if err != nil {
		fmt.Printf("Error: %v\n", err)
		return
	}
	fmt.Printf("%v", coefs.Result)
	// Output: [2 1 2 2 1 2]
}
