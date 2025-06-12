package gosynthcalc

import (
	"fmt"
	"io"
	"log"
	"os"
	"strings"
	"testing"
	"time"
)

func setup(fname string) ([]ChemicalReaction, []ChemicalFormula) {
	file, err := os.Open(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer func() {
		if err = file.Close(); err != nil {
			log.Fatal(err)
		}
	}()
	b, err := io.ReadAll(file)
	if err != nil {
		log.Fatal(err)
	}
	reactionsStr := strings.Split(string(b), "\n")
	reactions := []ChemicalReaction{}
	formulas := []ChemicalFormula{}
	for _, reac := range reactionsStr {
		if strings.TrimSpace(reac) == "" {
			continue
		}
		reacO, err := NewChemicalReaction(reac)
		if err != nil {
			panic(err)
		}
		reactions = append(reactions, *reacO)
		forms, err := reacO.ChemFormulas()
		if err != nil {
			panic(err)
		}
		formulas = append(formulas, forms...)
	}
	return reactions, formulas
}

func BenchmarkChemicalFormula_output(b *testing.B) {
	_, formulas := setup("data/text_mined_reactions.txt")

	f, err := os.Create("data/formula_output.txt")
	if err != nil {
		b.Fatal(err)
	}
	defer f.Close()

	var totalTime time.Duration
	var calls int64

	for b.Loop() {
		for _, form := range formulas {
			start := time.Now()
			out := form.Output()
			_, err := f.WriteString(out.String() + "\n")
			totalTime += time.Since(start)
			calls++
			if err != nil {
				b.Fatal(err)
			}
		}
	}

	if calls > 0 {
		avgTimeMilli := float64(totalTime.Milliseconds()) / float64(calls)
		b.ReportMetric(avgTimeMilli, "ms/formula")
	}
}

func BenchmarkChemicalReaction_output(b *testing.B) {
	reactions, _ := setup("data/text_mined_reactions.txt")

	f, err := os.Create("data/reaction_output.txt")
	if err != nil {
		b.Fatal(err)
	}
	defer f.Close()

	var totalTime time.Duration
	var calls int64

	for b.Loop() {
		for idx, reac := range reactions {
			start := time.Now()
			out, oerr := reac.Output()
			if oerr != nil {
				b.Fatalf("Error at reaction %d: %v", idx, oerr)
			}
			_, err := f.WriteString(fmt.Sprintf("%s\n", out))
			totalTime += time.Since(start)
			calls++
			if err != nil {
				b.Fatal(err)
			}
		}
	}

	if calls > 0 {
		avgTimeMilli := float64(totalTime.Milliseconds()) / float64(calls)
		b.ReportMetric(avgTimeMilli, "ms/reaction")
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
