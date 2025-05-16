package chemreaction

import (
	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"gonum.org/v1/gonum/mat"
)

func createReacMatrix(parsedFormulas [][]chemformula.Atom) *mat.Dense {
	matrix := mat.NewDense(3, 5, nil)
	return matrix
}
