package chemformula

import (
	"fmt"
	"slices"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

type oxide struct {
	metal   string
	formula string
	massP   float64
}

type molarMass struct {
	parsed []Atom
}

func (m molarMass) atomicMasses() []float64 {
	masses := make([]float64, len(m.parsed))
	for i, atom := range m.parsed {
		masses[i] = periodicTable[atom.Label].weight * atom.Amount
	}
	return masses
}

func (m molarMass) molarMass() float64 {
	return utils.SumFloatS(m.atomicMasses())
}

func (m molarMass) massPercent() []Atom {
	percent := make([]Atom, len(m.parsed))
	atomicMasses := m.atomicMasses()
	molarMass := m.molarMass()
	for i, mass := range atomicMasses {
		percent[i] = Atom{Label: m.parsed[i].Label, Amount: mass / molarMass * 100}
	}
	return percent
}

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

func (m molarMass) customOxides(inOxides ...string) ([]oxide, error) {
	oxides := []oxide{}
	metals := []string{}
	for _, cOxide := range inOxides {
		validator := formulaValidator{formula: cOxide}
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
				label = periodicTable[atom.Label].defaultOxide
			}
			oxides = append(oxides, oxide{metal: atom.Label, formula: label, massP: massPercents[i].Amount})
		}

	}

	return oxides, nil
}

func (m molarMass) oxidePercent(inOxides ...string) ([]Atom, error) {
	ret := []Atom{}
	oxides, err := m.customOxides(inOxides...)
	if err != nil {
		return nil, err
	}

	oxPercents := []float64{}
	for _, oxide := range oxides {
		parsedOxide := ChemicalFormulaParser{}.parse(oxide.formula)
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
