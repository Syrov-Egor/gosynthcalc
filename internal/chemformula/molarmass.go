package chemformula

import (
	"fmt"
	"slices"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

// oxide describes one non-oxygen element of a compound together with the
// oxide its mass percent is expressed as: metal is the element symbol
// ("Fe"), formula is the oxide formula used as its label ("Fe2O3") and massP
// is the mass percent of that element in the original compound.
type oxide struct {
	metal   string
	formula string
	massP   float64
}

// molarMass computes the molar mass and the composition percentages of a
// parsed formula.
type molarMass struct {
	parsed []Atom
}

// atomicMasses returns the standard atomic weight of each element multiplied
// by its amount in the formula, i.e. each element's contribution to the
// molar mass.
func (m molarMass) atomicMasses() []float64 {
	masses := make([]float64, len(m.parsed))
	for i, atom := range m.parsed {
		masses[i] = periodicTable[atom.Label].weight * atom.Amount
	}
	return masses
}

// molarMass returns the sum of the atomic masses, i.e. the
// [molar mass](https://en.wikipedia.org/wiki/Molar_mass) of the compound in
// g/mol.
func (m molarMass) molarMass() float64 {
	return utils.SumFloatS(m.atomicMasses())
}

// massPercent returns the
// [mass percent](https://en.wikipedia.org/wiki/Mass_fraction_(chemistry)) of
// every element: its share of the molar mass in percent, with the values of
// all elements summing to 100.
func (m molarMass) massPercent() []Atom {
	percent := make([]Atom, len(m.parsed))
	atomicMasses := m.atomicMasses()
	molarMass := m.molarMass()
	for i, mass := range atomicMasses {
		percent[i] = Atom{Label: m.parsed[i].Label, Amount: mass / molarMass * 100}
	}
	return percent
}

// atomicPercent returns the
// [atomic percent](https://en.wikipedia.org/wiki/Mole_fraction) of every
// element: its share of the total amount of substance in percent, with the
// values of all elements summing to 100.
func (m molarMass) atomicPercent() []Atom {
	percent := make([]Atom, len(m.parsed))
	amounts := []float64{}
	for _, atom := range m.parsed {
		amounts = append(amounts, atom.Amount)
	}
	sum := utils.SumFloatS(amounts)
	for i, amount := range amounts {
		percent[i] = Atom{Label: m.parsed[i].Label, Amount: amount / sum * 100}
	}
	return percent
}

// customOxides pairs every non-oxygen element of the compound with the oxide
// formula its mass percent should be reported as: an oxide listed in
// inOxides wins over the element's default oxide from the periodic table.
//
// Each custom oxide must be a binary compound whose second element is
// oxygen, otherwise an error is returned.
func (m molarMass) customOxides(inOxides ...string) ([]oxide, error) {
	oxides := []oxide{}
	metals := []string{}
	for _, cOxide := range inOxides {
		formula, err := NewChemicalFormula(cOxide)
		if err != nil {
			return nil, err
		}

		parsed := formula.ParsedFormula()
		if len(parsed) != 2 {
			return nil, fmt.Errorf("Only binary compounds can be considered as input (oxide '%s')", cOxide)
		} else if parsed[1].Label != "O" {
			return nil, fmt.Errorf("Only oxides can be considered as input (oxide '%s')", cOxide)
		}

		metals = append(metals, parsed[0].Label)
	}

	cOxides := make(map[string]string)
	for i := range metals {
		cOxides[metals[i]] = inOxides[i]
	}

	massPercents := m.massPercent()
	label := ""
	for i, atom := range m.parsed {
		if atom.Label != "O" {
			if slices.Contains(metals, atom.Label) {
				label = cOxides[atom.Label]
			} else {
				label = periodicTable[atom.Label].defaultOxide
			}
			oxides = append(oxides, oxide{metal: atom.Label, formula: label, massP: massPercents[i].Amount})
		}

	}

	return oxides, nil
}

// oxidePercent returns the percentage of every element expressed as its
// oxide (oxygen itself is skipped), normalized to 100. This kind of data is
// mostly used in XRF spectrometry and mineralogy. Each element's mass
// percent is converted to its oxide percent with the
// [conversion factor](https://www.geol.umd.edu/~piccoli/probe/molweight.html)
// between the element and its oxide (molar mass of the oxide divided by the
// atomic weight of the element and by the element's coefficient in the
// oxide), and the results are then scaled so that they sum to 100.
//
// inOxides may contain custom oxide formulas that replace the defaults of
// the periodic table.
func (m molarMass) oxidePercent(inOxides ...string) ([]Atom, error) {
	ret := []Atom{}
	oxides, err := m.customOxides(inOxides...)
	if err != nil {
		return nil, err
	}

	oxPercents := []float64{}
	for _, oxide := range oxides {
		formula, err := NewChemicalFormula(oxide.formula)
		if err != nil {
			return nil, err
		}
		parsedOxide := formula.ParsedFormula()
		oxideMass := molarMass{parsedOxide}.molarMass()
		atomicOxideCoef := parsedOxide[0].Amount
		atomicMass := periodicTable[oxide.metal].weight
		convFactor := oxideMass / atomicMass / atomicOxideCoef
		oxPercents = append(oxPercents, oxide.massP*convFactor)
	}

	normOxPercents := []float64{}
	sumOxPercents := utils.SumFloatS(oxPercents)
	for _, percent := range oxPercents {
		normOxPercents = append(normOxPercents, percent/sumOxPercents*100)
	}
	for i, oxide := range oxides {
		ret = append(ret, Atom{Label: oxide.formula, Amount: normOxPercents[i]})
	}

	return ret, nil
}
