package chemreaction

import (
	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"gonum.org/v1/gonum/mat"
)

func createReacMatrix(parsedFormulas [][]chemformula.Atom) *mat.Dense {
	atomMap := make(map[string]int)
	var atomOrder []string

	for _, formula := range parsedFormulas {
		for _, atom := range formula {
			if _, exists := atomMap[atom.Label]; !exists {
				atomOrder = append(atomOrder, atom.Label)
				atomMap[atom.Label] = len(atomOrder) - 1
			}
		}
	}

	numAtoms := len(atomOrder)
	numFormulas := len(parsedFormulas)
	data := make([]float64, numAtoms*numFormulas)
	for formulaIdx, formula := range parsedFormulas {
		for _, atom := range formula {
			atomIdx := atomMap[atom.Label]
			dataIdx := atomIdx*numFormulas + formulaIdx
			data[dataIdx] = atom.Amount
		}
	}

	return mat.NewDense(numAtoms, numFormulas, data)
}
