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

	parsedFormula []Atom
	molarMass     *float64
	massPercent   *[]Atom
	atomicPercent *[]Atom
	oxidePercent  map[string][]Atom
}

func NewChemicalFormula(formula string, precision ...uint) (*ChemicalFormula, error) {
	var prec uint = DefaultPrecision
	if len(precision) > 0 {
		prec = precision[0]
	}

	if err := validateBrackets(formula); err != nil {
		return nil, err
	}
	sanFormula := sanitize(formula)
	parsed, err := NewParser(sanFormula).parse()
	if err != nil {
		return nil, err
	}

	return &ChemicalFormula{
		formula:       formula,
		sanFormula:    sanFormula,
		precision:     prec,
		parsedFormula: parsed,
	}, nil
}

func (c *ChemicalFormula) Formula() string {
	return c.formula
}

func (c *ChemicalFormula) ParsedFormula() []Atom {
	return c.parsedFormula
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
	key := oxideCacheKey(inOxides)
	if cached, ok := c.oxidePercent[key]; ok {
		return cached, nil
	}

	percent, err := molarMass{c.ParsedFormula()}.oxidePercent(inOxides...)
	if err != nil {
		return nil, err
	}
	percent = roundAtomS(percent, c.precision)
	if c.oxidePercent == nil {
		c.oxidePercent = make(map[string][]Atom)
	}
	c.oxidePercent[key] = percent
	return percent, nil
}

func oxideCacheKey(inOxides []string) string {
	return strings.Join(inOxides, "\x00")
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
		Formula:        c.formula,
		ParsedFormula:  c.ParsedFormula(),
		MolarMass:      c.MolarMass(),
		MassPercent:    c.MassPercent(),
		AtomicPercent:  c.AtomicPercent(),
		OxidePercent:   oxides,
		printPrecision: pPrecision,
	}

	return cfO
}

type cfOutput struct {
	Formula        string
	ParsedFormula  []Atom
	MolarMass      float64
	MassPercent    []Atom
	AtomicPercent  []Atom
	OxidePercent   []Atom
	printPrecision uint
}

func (o cfOutput) String() string {
	var sb strings.Builder
	fmt.Fprintln(&sb, "formula:", o.Formula)
	fmt.Fprintln(&sb, "parsed formula:", o.ParsedFormula)
	fmt.Fprintf(&sb, "molar mass: %.*f\n", o.printPrecision, o.MolarMass)
	fmt.Fprintln(&sb, "mass percent:", formatAtoms(o.MassPercent, o.printPrecision))
	fmt.Fprintln(&sb, "atomic percent:", formatAtoms(o.AtomicPercent, o.printPrecision))
	fmt.Fprint(&sb, "oxide percent: ", formatAtoms(o.OxidePercent, o.printPrecision))
	return sb.String()
}

func formatAtoms(atoms []Atom, precision uint) string {
	var sb strings.Builder
	sb.WriteByte('[')
	for i, atom := range atoms {
		if i > 0 {
			sb.WriteByte(' ')
		}
		fmt.Fprintf(&sb, "'%s': %.*f", atom.Label, precision, atom.Amount)
	}
	sb.WriteByte(']')
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
