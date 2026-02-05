package chemformula

import (
	"unicode"
)

type formulaValidator struct {
	initialFormula string
	sanFormula     string
}

type validationResults struct {
	invalidSymbols []rune
	invalidAtoms   []string
	lBracket       int
	rBracket       int
	adduct         int
	nOfLetters     int
}

func (v formulaValidator) validate() error {
	if len(v.sanFormula) == 0 {
		return EmptyFormulaError{}
	}

	r := v.validatorRunner()

	if r.nOfLetters == 0 {
		return NoLettersError{v.initialFormula}
	}

	if len(r.invalidSymbols) > 0 {
		return InvalidSymbolsError{v.initialFormula, r.invalidSymbols}
	}

	if len(r.invalidAtoms) > 0 {
		return InvalidAtomsError{v.initialFormula, r.invalidAtoms}
	}

	if r.lBracket != r.rBracket {
		return BracketsNotBalancedError{v.initialFormula}
	}

	if r.adduct > 1 {
		return MoreThanOneAdductError{v.initialFormula}
	}

	return nil
}

func (v formulaValidator) validatorRunner() validationResults {
	invalidSymbols := make([]rune, 0, len(v.sanFormula))
	invalidAtoms := make([]string, 0, len(v.sanFormula))
	var lBracket, rBracket, adduct, letters int

	for i := 0; i < len(v.sanFormula); i++ {
		char := rune(v.sanFormula[i])

		if !v.isSymbolAllowed(char) {
			invalidSymbols = append(invalidSymbols, char)
		}

		switch char {
		case '(':
			lBracket++
		case ')':
			rBracket++
		case '*':
			adduct++
		}

		if unicode.IsLetter(char) {
			letters++
		}

		if unicode.IsLower(char) {
			if i == 0 || !unicode.IsUpper(rune(v.sanFormula[i-1])) {
				invalidAtoms = append(invalidAtoms, string(char))
			}
		}

		if unicode.IsUpper(char) {
			elem := string(v.sanFormula[i])

			if i+1 < len(v.sanFormula) && unicode.IsLower(rune(v.sanFormula[i+1])) {
				elem = v.sanFormula[i : i+2]
				i++
			}

			_, ok := periodicTable[elem]

			if !ok {
				invalidAtoms = append(invalidAtoms, elem)
			}
		}
	}
	return validationResults{
		invalidSymbols: invalidSymbols,
		invalidAtoms:   invalidAtoms,
		lBracket:       lBracket,
		rBracket:       rBracket,
		adduct:         adduct,
		nOfLetters:     letters,
	}
}

func (v formulaValidator) isSymbolAllowed(r rune) bool {
	return (r >= 'A' && r <= 'Z') ||
		(r >= 'a' && r <= 'z') ||
		(r >= '0' && r <= '9') ||
		r == '(' || r == ')' ||
		r == '*' || r == '.'
}
