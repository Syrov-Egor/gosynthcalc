package chemformula

import (
	"fmt"
	"slices"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

type Oxide struct {
	metal   string
	formula string
	massP   float64
}

type MolarMass struct {
	parsed []Atom
}

func (m MolarMass) atomicMasses() []float64 {
	masses := make([]float64, len(m.parsed))
	for i, atom := range m.parsed {
		masses[i] = PeriodicTable[atom.Label].Weight * atom.Amount
	}
	return masses
}

func (m MolarMass) molarMass() float64 {
	return utils.SumFloatS(m.atomicMasses())
}

func (m MolarMass) massPercent() []Atom {
	percent := make([]Atom, len(m.parsed))
	atomicMasses := m.atomicMasses()
	molarMass := m.molarMass()
	for i, mass := range atomicMasses {
		percent[i] = Atom{Label: m.parsed[i].Label, Amount: mass / molarMass * 100}
	}
	return percent
}

func (m MolarMass) atomicPercent() []Atom {
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

func (m MolarMass) customOxides(inOxides ...string) ([]Oxide, error) {
	oxides := []Oxide{}
	metals := []string{}
	for _, cOxide := range inOxides {
		validator := FormulaValidator{formula: cOxide}
		err := validator.validate()
		if err != nil {
			return nil, err
		}

		parsed := ChemicalFormulaParser{}.parse(cOxide)
		if len(parsed) > 2 {
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
				label = PeriodicTable[atom.Label].DefaultOxide
			}
			oxides = append(oxides, Oxide{metal: atom.Label, formula: label, massP: massPercents[i].Amount})
		}

	}

	return oxides, nil
}

func (m MolarMass) oxidePercent(inOxides ...string) ([]Atom, error) {
	ret := []Atom{}
	oxides, err := m.customOxides(inOxides...)
	if err != nil {
		return nil, err
	}

	oxPercents := []float64{}
	for _, oxide := range oxides {
		parsedOxide := ChemicalFormulaParser{}.parse(oxide.formula)
		oxideMass := MolarMass{parsedOxide}.molarMass()
		atomicOxideCoef := parsedOxide[0].Amount
		atomicMass := PeriodicTable[oxide.metal].Weight
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
