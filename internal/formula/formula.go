package formula

import "strings"

type ChemicalFormula struct {
	formula       string
	parsedFormula *[]Atom
	molarMass     *float64
}

func NewChemicalFormula(formula string) (*ChemicalFormula, error) {
	validator := FormulaValidator{formula: formula}
	err := validator.validate()
	if err != nil {
		return nil, err
	}
	newFormula := strings.Replace(formula, " ", "", -1)
	return &ChemicalFormula{
		formula: newFormula,
	}, nil
}

func (c *ChemicalFormula) Formula() string {
	return c.formula
}

func (c *ChemicalFormula) ParsedFormula() []Atom {
	if c.parsedFormula == nil {
		parser := ChemicalFormulaParser{}
		parsed := parser.parse(c.formula)
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
