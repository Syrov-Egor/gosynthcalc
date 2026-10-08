// Package chemformula parses and analyzes single chemical formulas.
//
// The entry point is [ChemicalFormula]: a formula string is validated and
// parsed once, and the parsed result is then used to
// compute the molar mass, mass percent, atomic percent and oxide percent of
// the compound.
//
// # Accepted input
//
// A formula may contain element symbols (canonical case, e.g. "Fe"), digits
// and decimal points for fractional atom amounts, nested parentheses in any
// of the three styles "()", "[]" and "{}" (all treated alike), and a single
// adduct separator written as "*", "·" or "•", as in "CuSO4*5H2O".
// Whitespace is ignored. Any other character, unknown element symbol,
// unbalanced bracket or second adduct symbol makes construction fail with an
// error describing the problem (the parser additionally reports the 1-based
// position of the offending token within the normalized formula).
//
// All calculated values are rounded to the precision of the formula
// ([DefaultPrecision] unless overridden).
package chemformula

import (
	"fmt"
	"strings"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

const (
	// DefaultPrecision is the number of decimal places used to round
	// calculated values of a [ChemicalFormula].
	DefaultPrecision uint = 8
	// DefaultPrintPrecision is the number of decimal places used when
	// rendering results with [ChemicalFormula.Output].
	DefaultPrintPrecision uint = 4
)

// ChemicalFormula represents a single chemical formula and lazily computes
// properties from it. Construct it with [NewChemicalFormula], which
// validates and parses the string; every property below is then calculated
// on first use and cached for the lifetime of the value.
//
// The slices returned by the accessors are the cached values themselves:
// callers must not modify them.
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

// NewChemicalFormula validates formula, normalizes it and parses it into
// atoms. The optional precision sets the rounding precision of every
// calculated property and defaults to [DefaultPrecision].
//
// It fails if the brackets are unbalanced or mismatched, if the formula is
// empty, contains characters outside the accepted set or unknown element
// symbols, contains malformed numbers or more than one adduct symbol.
//
// Example:
//
//	f, err := NewChemicalFormula("H2O")
//	f.MolarMass()   // 18.015
//	f.MassPercent() // ['H': 11.19067444 'O': 88.80932556]
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

// Formula returns the formula string as it was passed to
// [NewChemicalFormula], whitespace included; normalization happens only
// internally.
func (c *ChemicalFormula) Formula() string {
	return c.formula
}

// ParsedFormula returns the formula as atoms ordered by the first appearance
// of each element in the string.
func (c *ChemicalFormula) ParsedFormula() []Atom {
	return c.parsedFormula
}

// MolarMass returns the
// [molar mass](https://en.wikipedia.org/wiki/Molar_mass) of the compound in
// g/mol, rounded to the formula's precision.
func (c *ChemicalFormula) MolarMass() float64 {
	if c.molarMass == nil {
		mass := molarMass{c.ParsedFormula()}.molarMass()
		mass = utils.RoundFloat(mass, c.precision)
		c.molarMass = &mass
	}
	return *c.molarMass
}

// MassPercent returns the
// [mass percent](https://en.wikipedia.org/wiki/Mass_fraction_(chemistry)) of
// every element of the compound: values in percent that sum to 100, rounded
// to the formula's precision.
func (c *ChemicalFormula) MassPercent() []Atom {
	if c.massPercent == nil {
		percent := molarMass{c.ParsedFormula()}.massPercent()
		percent = roundAtomS(percent, c.precision)
		c.massPercent = &percent
	}
	return *c.massPercent
}

// AtomicPercent returns the
// [atomic percent](https://en.wikipedia.org/wiki/Mole_fraction) of every
// element of the compound: mole fractions in percent that sum to 100,
// rounded to the formula's precision.
func (c *ChemicalFormula) AtomicPercent() []Atom {
	if c.atomicPercent == nil {
		percent := molarMass{c.ParsedFormula()}.atomicPercent()
		percent = roundAtomS(percent, c.precision)
		c.atomicPercent = &percent
	}
	return *c.atomicPercent
}

// OxidePercent returns the percentage of every oxide in the compound,
// normalized to 100: each element except oxygen is reported as the oxide it
// normally forms, whose formula comes from the periodic table unless a
// custom one is supplied. Custom oxides may be passed on each call; they must
// be binary compounds ending with oxygen, otherwise an error is returned.
// Results are cached per set of custom oxides.
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

// oxideCacheKey builds the cache key of an oxide percent computation out of
// the custom oxide formulas (NUL-separated, so no formula can alias another).
func oxideCacheKey(inOxides []string) string {
	return strings.Join(inOxides, "\x00")
}

// Output gathers every property of the formula into a [cfOutput] value in a
// single pass; it corresponds to the output_results dict of the Python class.
// The optional printPrecision overrides [DefaultPrintPrecision] and controls
// only how many digits [cfOutput.String] prints, not the stored values.
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

// cfOutput is the collected result of a [ChemicalFormula] computation, the Go
// counterpart of the output_results dict of the Python class. Its
// [cfOutput.String] method renders the human-readable report.
type cfOutput struct {
	Formula        string
	ParsedFormula  []Atom
	MolarMass      float64
	MassPercent    []Atom
	AtomicPercent  []Atom
	OxidePercent   []Atom
	printPrecision uint
}

// String formats the report: one line per property, with the percentage tables
// and the molar mass rounded to printPrecision digits.
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

// formatAtoms renders a percent table as "'K': 44.8752 'S': 18.3986 ..."
// with every value rounded to precision.
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

// roundAtomS returns a copy of s with every amount rounded to precision.
func roundAtomS(s []Atom, precision uint) []Atom {
	res := make([]Atom, len(s))
	for i, atom := range s {
		res[i] = Atom{Label: atom.Label,
			Amount: utils.RoundFloat(atom.Amount, precision)}
	}
	return res
}
