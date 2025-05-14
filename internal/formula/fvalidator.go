package formula

import (
	"fmt"
	"slices"
	"sort"
	"strings"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

type FormulaValidator struct {
	formula string
}

func (v FormulaValidator) emptyFormula() bool {
	return v.formula == ""
}

func (v FormulaValidator) invalidCharacters() []string {
	return regexes.allowedSymbols.FindAllString(v.formula, -1)
}

func (v FormulaValidator) invalidAtoms() []string {
	atoms := regexes.atomRegex.FindAllString(v.formula, -1)
	invalid := make([]string, 0)
	cFormula := strings.Clone(v.formula)
	slices.Sort(atoms)
	uniqueAtoms := slices.Compact(atoms)
	sort.Slice(uniqueAtoms, func(i, j int) bool {
		return len(uniqueAtoms[i]) > len(uniqueAtoms[j])
	})

	for _, atom := range uniqueAtoms {
		if !slices.Contains(PeriodicTableElements, atom) {
			invalid = append(invalid, atom)
		}
		cFormula = strings.Replace(cFormula, atom, "", -1)
	}
	leftovers := regexes.letterRegex.FindAllString(cFormula, -1)
	invalid = append(invalid, leftovers...)
	return invalid
}

func (v FormulaValidator) bracketsBalance() bool {
	counter := utils.StringCounter(v.formula)
	for i := range len(regexes.openerBrackets) {
		open := string(regexes.openerBrackets[i])
		close := string(regexes.closerBrackets[i])
		if counter[open] != counter[close] {
			return false
		}
	}
	return true
}

func (v FormulaValidator) numOfAdducts() int {
	counter := utils.StringCounter(v.formula)
	i := 0
	for _, adduct := range regexes.adductSymbols {
		i += counter[string(adduct)]
	}
	return i
}

func (v FormulaValidator) validate() error {
	var err error
	switch {
	case v.emptyFormula():
		err = fmt.Errorf("Empty formula string")
	case len(v.invalidCharacters()) > 0:
		err = fmt.Errorf("There are invalid characters in the formula: %s", v.invalidCharacters())
	case len(v.invalidAtoms()) > 0:
		err = fmt.Errorf("There are invalid atoms in the formula: %s", v.invalidAtoms())
	case !v.bracketsBalance():
		err = fmt.Errorf("Brackets %s %s are not balanced", string(regexes.openerBrackets), string(regexes.closerBrackets))
	case v.numOfAdducts() > 1:
		err = fmt.Errorf("There are more than 1 adduct symbol %s", string(regexes.adductSymbols))
	}
	return err
}
