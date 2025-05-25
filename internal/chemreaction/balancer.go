package chemreaction

import "gonum.org/v1/gonum/mat"

type Balancer struct {
	reactionMatrix *mat.Dense
	separatorPos   int
	intify         bool
	tolerance      float64
	bAlgos         *BalancingAlgos
}

func NewBalancer(matrix *mat.Dense, separatorPos int, intify bool, tolerance ...float64) *Balancer {
	var tol float64
	if tolerance == nil {
		tol = 1e-12
	} else {
		tol = tolerance[0]
	}

	bAlgos := NewBalancingAlgos(matrix, separatorPos, tol)

	return &Balancer{
		reactionMatrix: matrix,
		separatorPos:   separatorPos,
		intify:         intify,
		tolerance:      tol,
		bAlgos:         bAlgos,
	}
}

func (b *Balancer) intifyCoefs(coefs []float64, limit int) {

}
