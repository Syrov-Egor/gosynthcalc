package chemformula

import (
	"fmt"
	"strings"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

const (
	DefaultPrecision      uint = 8
	DefaultPrintPrecision uint = 4
)

type ChemicalFormula struct {
	formula    string
	sanFormula string
	precision  uint

	parsedFormula *[]Atom
	molarMass     *float64
	massPercent   *[]Atom
	atomicPercent *[]Atom
	oxidePercent  *[]Atom
}

func NewChemicalFormula(formula string, precision ...uint) (*ChemicalFormula, error) {
	var prec uint = DefaultPrecision
	if len(precision) > 0 {
		prec = precision[0]
	}

	sanFormula := sanitize(formula)
	validator := formulaValidator{formula: sanFormula}
	err := validator.validate()
	if err != nil {
		return nil, err
	}

	return &ChemicalFormula{
		formula:    formula,
		sanFormula: sanFormula,
		precision:  prec,
	}, nil
}

func (c *ChemicalFormula) Formula() string {
	return c.formula
}

func (c *ChemicalFormula) ParsedFormula() []Atom {
	if c.parsedFormula == nil {
		parser := NewParser(c.sanFormula)
		parsed := parser.parse()
		c.parsedFormula = &parsed
	}
	return *c.parsedFormula
}

func (c *ChemicalFormula) MolarMass() float64 {
	if c.molarMass == nil {
		mass := molarMass{c.ParsedFormula()}.molarMass()
		mass = utils.RoundFloat(mass, c.precision)
		c.molarMass = &mass
	}
	return *c.molarMass
}

func (c *ChemicalFormula) MassPercent() []Atom {
	if c.massPercent == nil {
		percent := molarMass{c.ParsedFormula()}.massPercent()
		percent = roundAtomS(percent, c.precision)
		c.massPercent = &percent
	}
	return *c.massPercent
}

func (c *ChemicalFormula) AtomicPercent() []Atom {
	if c.atomicPercent == nil {
		percent := molarMass{c.ParsedFormula()}.atomicPercent()
		percent = roundAtomS(percent, c.precision)
		c.atomicPercent = &percent
	}
	return *c.atomicPercent
}

func (c *ChemicalFormula) OxidePercent(inOxides ...string) ([]Atom, error) {
	if c.oxidePercent == nil {
		percent, err := molarMass{c.ParsedFormula()}.oxidePercent(inOxides...)
		if err != nil {
			return nil, err
		}
		percent = roundAtomS(percent, c.precision)
		c.oxidePercent = &percent
	}
	return *c.oxidePercent, nil
}

func (c *ChemicalFormula) Output(printPrecision ...uint) cfOutput {
	var pPrecision uint
	if printPrecision == nil {
		pPrecision = DefaultPrintPrecision
	} else {
		pPrecision = printPrecision[0]
	}

	oxides, _ := c.OxidePercent()
	cfO := cfOutput{
		Formula:       c.formula,
		ParsedFormula: c.ParsedFormula(),
		MolarMass:     utils.RoundFloat(c.MolarMass(), pPrecision),
		MassPercent:   roundAtomS(c.MassPercent(), pPrecision),
		AtomicPercent: roundAtomS(c.AtomicPercent(), pPrecision),
		OxidePercent:  roundAtomS(oxides, pPrecision),
	}

	return cfO
}

type cfOutput struct {
	Formula       string
	ParsedFormula []Atom
	MolarMass     float64
	MassPercent   []Atom
	AtomicPercent []Atom
	OxidePercent  []Atom
}

func (o cfOutput) String() string {
	var sb strings.Builder
	fmt.Fprintln(&sb, "formula:", o.Formula)
	fmt.Fprintln(&sb, "parsed formula:", o.ParsedFormula)
	fmt.Fprintln(&sb, "molar mass:", o.MolarMass)
	fmt.Fprintln(&sb, "mass percent:", o.MassPercent)
	fmt.Fprintln(&sb, "atomic percent:", o.AtomicPercent)
	fmt.Fprint(&sb, "oxide percent: ", o.OxidePercent)
	return sb.String()
}

func roundAtomS(s []Atom, precision uint) []Atom {
	res := make([]Atom, len(s))
	for i, atom := range s {
		res[i] = Atom{Label: atom.Label,
			Amount: utils.RoundFloat(atom.Amount, precision)}
	}
	return res
}
