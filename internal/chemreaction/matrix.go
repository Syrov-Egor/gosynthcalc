package chemreaction

import (
	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"gonum.org/v1/gonum/mat"
)

// createReacMatrix builds the reaction matrix out of the parsed formulas:
// rows are the distinct elements in the order of their first appearance in
// the reaction, columns are the compounds, and the entry (i, j) holds the
// amount of element i in compound j. A zero means the element is absent from
// that compound.
//
// The first implementation of the reaction matrix method probably belongs to
// [Blakley](https://doi.org/10.1021/ed059p728). In general, a chemical
// reaction matrix is composed of the coefficients of each atom in each
// compound, giving a 2D array.
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
