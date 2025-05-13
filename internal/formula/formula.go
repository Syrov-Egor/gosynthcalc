package formula

type ChemicalFormula struct {
	formula       string
	parsedFormula *[]Atom
}

func NewChemicalFormula(formula string) *ChemicalFormula {
	return &ChemicalFormula{
		formula: formula,
	}
}

func (c *ChemicalFormula) ParsedFormula() []Atom {
	if c.parsedFormula == nil {
		parsed := NewChemicalFormulaParser().parse(c.formula)
		c.parsedFormula = &parsed
	}
	return *c.parsedFormula
}
