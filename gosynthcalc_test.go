package gosynthcalc

import (
	"fmt"
	"io"
	"log"
	"os"
	"strings"
	"testing"
)

func setup(fname string) ([]string, []string) {
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
	formulas := []string{}
	reactions := []string{}

	for _, reac := range reactionsStr {
		if strings.TrimSpace(reac) == "" {
			continue
		}
		reacO, err := NewChemicalReaction(reac)
		if err != nil {
			panic(err)
		}
		reactions = append(reactions, reac)

		forms, err := reacO.ChemFormulas()
		if err != nil {
			panic(err)
		}
		for _, f := range forms {
			formulas = append(formulas, f.Formula())
		}
	}
	return formulas, reactions
}

func BenchmarkChemicalFormula_output(b *testing.B) {
	b.ReportAllocs()
	formulas, _ := setup("data/text_mined_reactions.txt")
	n := len(formulas)

	i := 0
	for b.Loop() {
		form := formulas[i%n]
		i++

		formulaObj, err := NewChemicalFormula(form)
		if err != nil {
			b.Fatal(err)
		}
		_ = formulaObj.Output()
	}
}

func BenchmarkChemicalReaction_output(b *testing.B) {
	b.ReportAllocs()
	_, reactions := setup("data/text_mined_reactions.txt")
	n := len(reactions)

	i := 0
	for b.Loop() {
		reac := reactions[i%n]
		i++

		reactionObj, err := NewChemicalReaction(reac)
		if err != nil {
			b.Fatal(err)
		}
		if _, err := reactionObj.Output(); err != nil {
			b.Fatal(err)
		}
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
