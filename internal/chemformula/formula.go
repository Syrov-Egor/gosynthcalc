package chemformula

import (
	"fmt"
	"strings"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

// TODO! "[]" formula testcase

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
	validator := formulaValidator{formula: newFormula}
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
		parser := chemicalFormulaParser{}
		parsed := parser.parse(c.formula)
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
		pPrecision = 4
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
	form := fmt.Sprintln("formula:", o.Formula)
	pForm := fmt.Sprintln("parsed formula:", o.ParsedFormula)
	mMass := fmt.Sprintln("molar mass:", o.MolarMass)
	mPercent := fmt.Sprintln("mass percent:", o.MassPercent)
	aPercent := fmt.Sprintln("atomic percent:", o.AtomicPercent)
	oPercent := fmt.Sprint("oxide percent: ", o.OxidePercent)
	res := form + pForm + mMass + mPercent + aPercent + oPercent
	return res
}

func roundAtomS(s []Atom, precision uint) []Atom {
	ret := make([]Atom, len(s))
	for i, atom := range s {
		ret[i] = Atom{Label: atom.Label,
			Amount: utils.RoundFloat(atom.Amount, precision)}
	}
	return ret
}
