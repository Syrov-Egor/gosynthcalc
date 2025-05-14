package formula

import (
	"fmt"
	"strings"
)

type ChemicalFormula struct {
	formula       string
	parsedFormula *[]Atom
	molarMass     *float64
	massPercent   *[]Atom
	atomicPercent *[]Atom
	oxidePercent  *[]Atom
}

func NewChemicalFormula(formula string) (*ChemicalFormula, error) {
	newFormula := strings.Replace(formula, " ", "", -1)
	validator := FormulaValidator{formula: newFormula}
	err := validator.validate()
	if err != nil {
		return nil, err
	}
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
		mass := MolarMass{c.ParsedFormula()}.molarMass()
		c.molarMass = &mass
	}
	return *c.molarMass
}

func (c *ChemicalFormula) MassPercent() []Atom {
	if c.massPercent == nil {
		percent := MolarMass{c.ParsedFormula()}.massPercent()
		c.massPercent = &percent
	}
	return *c.massPercent
}

func (c *ChemicalFormula) AtomicPercent() []Atom {
	if c.atomicPercent == nil {
		percent := MolarMass{c.ParsedFormula()}.atomicPercent()
		c.atomicPercent = &percent
	}
	return *c.atomicPercent
}

func (c *ChemicalFormula) OxidePercent() []Atom {
	if c.oxidePercent == nil {
		fmt.Println(MolarMass{c.ParsedFormula()}.customOxides("SO4", "K2O2"))
	}
	return []Atom{}
}
