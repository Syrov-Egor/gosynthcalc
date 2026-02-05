package chemformula

import "fmt"

type EmptyFormulaError struct {
}

func (e EmptyFormulaError) Error() string {
	return "Empty formula"
}

type NoLettersError struct {
	formula string
}

func (e NoLettersError) Error() string {
	return fmt.Sprintf("No letters A-Z or a-z in the formula '%s'",
		e.formula)
}

type InvalidSymbolsError struct {
	formula string
	symbols []rune
}

func (e InvalidSymbolsError) Error() string {
	invalid := make([]string, len(e.symbols))
	for i := range e.symbols {
		invalid[i] = string(e.symbols[i])
	}
	return fmt.Sprintf("There are invalid symbols(s) %v in the formula '%s'",
		invalid, e.formula)
}

type InvalidAtomsError struct {
	formula string
	atoms   []string
}

func (e InvalidAtomsError) Error() string {
	return fmt.Sprintf("There are invalid atoms(s) %v in the formula '%s'",
		e.atoms, e.formula)
}

type BracketsNotBalancedError struct {
	formula string
}

func (e BracketsNotBalancedError) Error() string {
	return fmt.Sprintf("Brackets ()[]{} are not balanced in the formula '%s'",
		e.formula)
}

type MoreThanOneAdductError struct {
	formula string
}

func (e MoreThanOneAdductError) Error() string {
	return fmt.Sprintf("There are more than 1 adduct symbol *·• in the formula '%s'",
		e.formula)
}
