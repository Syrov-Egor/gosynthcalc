package chemreaction

import (
	"fmt"

	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
)

type Coeffs struct {
	mode               Mode
	parsedFormulas     [][]chemformula.Atom
	decomposedReaction *reactionDecomposer
	balancer           *Balancer
}

func (c *Coeffs) calculateCoeffs() (MethodResult, error) {
	user := "user"
	switch c.mode {

	case force:
		return MethodResult{Method: user, Result: c.decomposedReaction.initCoefs}, nil

	case check:
		if isReactionBalanced(c.balancer.bAlgos.ReactantMatrix,
			c.balancer.bAlgos.ProductMatrix,
			c.decomposedReaction.initCoefs,
			c.balancer.tolerance) {
			return MethodResult{Method: user, Result: c.decomposedReaction.initCoefs}, nil
		} else {
			return MethodResult{Method: user, Result: nil}, fmt.Errorf("Reaction is not balanced")
		}

	case balance:
		coefs, err := c.balancer.Auto()
		if err != nil {
			return MethodResult{Method: user, Result: nil}, err
		}
		return coefs, nil

	default:
		return MethodResult{Method: user, Result: nil}, fmt.Errorf("No such mode %d", c.mode)
	}
}

func (c *Coeffs) validateCoeffs(coefs []float64) error {
	if !allPositive(coefs) {
		return fmt.Errorf("Some coefs in %v are negative", coefs)
	}
	_, cols := c.balancer.reactionMatrix.Dims()
	if len(coefs) != cols {
		return fmt.Errorf("Number of coefs should be equal %d, got %d", cols, len(coefs))
	}

	return nil
}
