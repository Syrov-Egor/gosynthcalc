package formula

type ChemicalFormula struct {
	formula       string
	parsedFormula *map[string]float64
}

func NewChemicalFormula(formula string) *ChemicalFormula {
	return &ChemicalFormula{
		formula: formula,
	}
}

func (c *ChemicalFormula) ParseFormula() map[string]float64 {
	if c.parsedFormula == nil {
		parsed, _ := NewChemicalFormulaParser().parse(c.formula)
		c.parsedFormula = &parsed
	}
	return *c.parsedFormula
}
