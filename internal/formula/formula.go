package formula

type ChemicalFormula struct {
	formula       string
	parsedFormula *[]Atom
	molarMass     *float64
	validation    *bool
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

func (c *ChemicalFormula) MolarMass() float64 {
	if c.molarMass == nil {
		mass := MolarMass{c.ParsedFormula()}.calcMolarMass()
		c.molarMass = &mass
	}
	return *c.molarMass
}

func (c *ChemicalFormula) Validation() bool {
	if c.validation == nil {
		validator := NewFormulaValidator(c.formula)
		validator.invalidAtoms()
	}
	return true
}
