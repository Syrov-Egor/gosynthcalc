package chemformula

import (
	"strings"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

type ChemicalFormula struct {
	formula       string
	precision     uint
	parsedFormula *[]Atom
	molarMass     *float64
	massPercent   *[]Atom
	atomicPercent *[]Atom
	oxidePercent  *[]Atom
}

func NewChemicalFormula(formula string, precision ...uint) (*ChemicalFormula, error) {
	var prec uint = 8
	if len(precision) > 0 {
		prec = precision[0]
	}

	newFormula := strings.Replace(formula, " ", "", -1)
	validator := FormulaValidator{formula: newFormula}
	err := validator.validate()
	if err != nil {
		return nil, err
	}

	return &ChemicalFormula{
		formula:   newFormula,
		precision: prec,
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
		mass = utils.RoundFloat(mass, c.precision)
		c.molarMass = &mass
	}
	return *c.molarMass
}

func (c *ChemicalFormula) MassPercent() []Atom {
	if c.massPercent == nil {
		percent := MolarMass{c.ParsedFormula()}.massPercent()
		percent = roundAtomS(percent, c.precision)
		c.massPercent = &percent
	}
	return *c.massPercent
}

func (c *ChemicalFormula) AtomicPercent() []Atom {
	if c.atomicPercent == nil {
		percent := MolarMass{c.ParsedFormula()}.atomicPercent()
		percent = roundAtomS(percent, c.precision)
		c.atomicPercent = &percent
	}
	return *c.atomicPercent
}

func (c *ChemicalFormula) OxidePercent(inOxides ...string) ([]Atom, error) {
	if c.oxidePercent == nil {
		percent, err := MolarMass{c.ParsedFormula()}.oxidePercent(inOxides...)
		if err != nil {
			return nil, err
		}
		percent = roundAtomS(percent, c.precision)
		c.oxidePercent = &percent
	}
	return *c.oxidePercent, nil
}

func roundAtomS(s []Atom, precision uint) []Atom {
	ret := make([]Atom, len(s))
	for i, atom := range s {
		ret[i] = Atom{Label: atom.Label,
			Amount: utils.RoundFloat(atom.Amount, precision)}
	}
	return ret
}
