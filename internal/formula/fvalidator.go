package formula

import (
	"fmt"
	"regexp"
	"slices"
	"sort"
	"strings"
)

type FormulaValidator struct {
	formula        string
	atomRegex      *regexp.Regexp
	allowedSymbols *regexp.Regexp
	letterRegex    *regexp.Regexp
}

func NewFormulaValidator(formula string) *FormulaValidator {
	return &FormulaValidator{
		formula:        formula,
		atomRegex:      regexp.MustCompile(`([A-Z][a-z]*)`),
		allowedSymbols: regexp.MustCompile(`[^A-Za-z0-9.({[)}\]*·•]`),
		letterRegex:    regexp.MustCompile(`[a-z]`),
	}
}

func (v FormulaValidator) emptyFormula() bool {
	return v.formula == ""
}

func (v FormulaValidator) invalidCharacters() []string {
	return v.allowedSymbols.FindAllString(v.formula, -1)
}

func (v FormulaValidator) invalidAtoms() []string {
	atoms := v.atomRegex.FindAllString(v.formula, -1)
	invalid := make([]string, 0)
	cFormula := strings.Clone(v.formula)

	seen := make(map[string]bool)
	uniqueAtoms := []string{}
	for _, atom := range atoms {
		if !seen[atom] {
			seen[atom] = true
			uniqueAtoms = append(uniqueAtoms, atom)
		}
	}
	sort.Slice(uniqueAtoms, func(i, j int) bool {
		return len(uniqueAtoms[i]) > len(uniqueAtoms[j])
	})

	for _, atom := range uniqueAtoms {
		if !slices.Contains(PeriodicTableElements, atom) {
			invalid = append(invalid, atom)
		}
		cFormula = strings.Replace(cFormula, atom, "", -1)
	}
	leftovers := v.letterRegex.FindAllString(cFormula, -1)
	invalid = append(invalid, leftovers...)
	fmt.Println(invalid)
	return invalid
}
