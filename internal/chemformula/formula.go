package chemformula

import "fmt"

type Atom struct {
	Label  string
	Amount float64
}

func (a Atom) String() string {
	return fmt.Sprintf("'%s': %v", a.Label, a.Amount)
}

type ChemicalFormula struct {
	initialFormula string
	sanFormula     string
	precision      int
}

func NewChemicalFormula(formula string, precision ...int) (*ChemicalFormula, error) {
	var prec int = 8
	if len(precision) > 0 && precision[0] > 0 {
		prec = precision[0]
	}

	sanFormula := formulaSanitizer{}.sanitize(formula)
	err := formulaValidator{formula, sanFormula}.validate()

	if err != nil {
		return nil, err
	}

	return &ChemicalFormula{
		initialFormula: formula,
		sanFormula:     sanFormula,
		precision:      prec,
	}, nil
}

func (c *ChemicalFormula) Formula() string {
	return c.initialFormula
}
