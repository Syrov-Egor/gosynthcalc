package chemformula

import (
	"fmt"
	"slices"
	"sort"
	"strings"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

type formulaValidator struct {
	formula string
}

func (v formulaValidator) emptyFormula() bool {
	return v.formula == ""
}

func (v formulaValidator) invalidCharacters() []string {
	return formRegexes.allowedSymbols.FindAllString(v.formula, -1)
}

func (v formulaValidator) invalidAtoms() []string {
	atoms := formRegexes.atomRegex.FindAllString(v.formula, -1)
	invalid := make([]string, 0)
	cFormula := strings.Clone(v.formula)
	uniqueAtoms := utils.UniqueElems(atoms)
	sort.Slice(uniqueAtoms, func(i, j int) bool {
		return len(uniqueAtoms[i]) > len(uniqueAtoms[j])
	})

	for _, atom := range uniqueAtoms {
		if !slices.Contains(periodicTableElements, atom) {
			invalid = append(invalid, atom)
		}
		cFormula = strings.Replace(cFormula, atom, "", -1)
	}
	leftovers := formRegexes.letterRegex.FindAllString(cFormula, -1)
	invalid = append(invalid, leftovers...)
	return invalid
}

func (v formulaValidator) bracketsBalance() bool {
	counter := utils.StringCounter(v.formula)
	for i := range len(formRegexes.openerBrackets) {
		open := string(formRegexes.openerBrackets[i])
		close := string(formRegexes.closerBrackets[i])
		if counter[open] != counter[close] {
			return false
		}
	}
	return true
}

func (v formulaValidator) numOfAdducts() int {
	counter := utils.StringCounter(v.formula)
	i := 0
	for _, adduct := range formRegexes.adductSymbols {
		i += counter[string(adduct)]
	}
	return i
}

func (v formulaValidator) validate() error {
	var err error
	switch {
	case v.emptyFormula():
		err = fmt.Errorf("Empty formula string")
	case len(v.invalidCharacters()) > 0:
		err = fmt.Errorf("There are invalid character(s) %s in the formula '%s'",
			v.invalidCharacters(), v.formula)
	case len(v.invalidAtoms()) > 0:
		err = fmt.Errorf("There are invalid atom(s) %s in the formula '%s'",
			v.invalidAtoms(), v.formula)
	case !v.bracketsBalance():
		err = fmt.Errorf("Brackets %s %s are not balanced in the formula '%s'",
			string(formRegexes.openerBrackets), string(formRegexes.closerBrackets), v.formula)
	case v.numOfAdducts() > 1:
		err = fmt.Errorf("There are more than 1 adduct symbol %s in the formula '%s'",
			string(formRegexes.adductSymbols), v.formula)
	}
	return err
}
