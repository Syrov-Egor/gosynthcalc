package main

import (
	"fmt"
	"io"
	"log"
	"os"
	"strings"
	"testing"
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
	reactionsStr := strings.Split(string(b), "\n")
	reactions := []ChemicalReaction{}
	formulas := []ChemicalFormula{}
	for _, reac := range reactionsStr {
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

var reactions, formulas = setup("data/text_mined_reactions.txt")

func BenchmarkChemincalFormula_output(b *testing.B) {
	for _, form := range formulas {
		out := form.Output()
		fmt.Println(out)
	}
}

func BenchmarkChemincalReaction_output(b *testing.B) {
	for i, reac := range reactions {
		_, oerr := reac.Output()
		if oerr != nil {
			fmt.Println(i, reac)
		}
	}
}
