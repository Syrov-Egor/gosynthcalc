package chemreaction

import (
	"fmt"
	"math/big"

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

// createExactReacMatrix builds the signed exact reaction matrix A = [R | -P]
// from exact parsed formulas. Rows and columns follow the same element and
// compound order as [createReacMatrix], so the exact matrix and the float
// matrix describe the same reaction. Reactant columns stay positive and
// product columns are negated.
func createExactReacMatrix(parsed [][]chemformula.ExactAtom, separatorPos int) (*exactMatrix, error) {
	numFormulas := len(parsed)
	if numFormulas == 0 {
		return nil, fmt.Errorf("%w: reaction has no compounds", ErrInvalidMatrix)
	}
	if separatorPos < 0 || separatorPos > numFormulas {
		return nil, fmt.Errorf("%w: separator position %d is outside [0,%d]",
			ErrInvalidMatrix, separatorPos, numFormulas)
	}

	atomMap := make(map[string]int)
	var atomOrder []string
	for _, formula := range parsed {
		for _, atom := range formula {
			if _, exists := atomMap[atom.Label]; !exists {
				atomOrder = append(atomOrder, atom.Label)
				atomMap[atom.Label] = len(atomOrder) - 1
			}
		}
	}
	if len(atomOrder) == 0 {
		return nil, fmt.Errorf("%w: reaction has no elements", ErrInvalidMatrix)
	}

	data := make([][]*big.Rat, len(atomOrder))
	for i := range data {
		data[i] = make([]*big.Rat, numFormulas)
		for j := range data[i] {
			data[i][j] = new(big.Rat)
		}
	}
	for formulaIdx, formula := range parsed {
		for _, atom := range formula {
			if atom.Amount == nil {
				return nil, fmt.Errorf("%w: nil exact amount for %q",
					ErrInvalidMatrix, atom.Label)
			}
			row := atomMap[atom.Label]
			if formulaIdx >= separatorPos {
				data[row][formulaIdx].Sub(data[row][formulaIdx], atom.Amount)
			} else {
				data[row][formulaIdx].Add(data[row][formulaIdx], atom.Amount)
			}
		}
	}

	return &exactMatrix{rows: len(atomOrder), cols: numFormulas, data: data}, nil
}
