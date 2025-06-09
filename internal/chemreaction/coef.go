package chemreaction

import (
	"fmt"

	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

type coeffs struct {
	mode               Mode
	parsedFormulas     [][]chemformula.Atom
	decomposedReaction *reactionDecomposer
	balancer           *balancer
}

func (c *coeffs) calculateCoeffs() (MethodResult, error) {
	user := "User"
	switch c.mode {

	case Force:
		return MethodResult{Method: user, Result: c.decomposedReaction.initCoefs}, nil

	case Check:
		if isReactionBalanced(c.balancer.bAlgos.ReactantMatrix,
			c.balancer.bAlgos.ProductMatrix,
			c.decomposedReaction.initCoefs,
			c.balancer.tolerance) {
			return MethodResult{Method: user, Result: c.decomposedReaction.initCoefs}, nil
		} else {
			return MethodResult{Method: user, Result: nil}, fmt.Errorf("Reaction is not balanced")
		}

	case Balance:
		coefs, err := c.balancer.Auto()
		if err != nil {
			return MethodResult{Method: user, Result: nil}, err
		}
		return coefs, nil

	default:
		return MethodResult{Method: user, Result: nil}, fmt.Errorf("No such mode %d", c.mode)
	}
}

func (c *coeffs) validateCoeffs(coefs []float64) error {
	_, cols := c.balancer.reactionMatrix.Dims()

	switch {
	case !allPositive(coefs):
		return fmt.Errorf("Some coefs in %v are negative or 0", coefs)
	case len(coefs) != cols:
		return fmt.Errorf("Number of coefs should be equal %d, got %d", cols, len(coefs))
	default:
		return nil
	}
}

func (c *coeffs) elementCountValidation() []string {
	if c.mode != Force {
		r := make([]string, 0)
		reactants := c.parsedFormulas[:c.balancer.separatorPos]
		for _, reac := range reactants {
			for _, atom := range reac {
				r = append(r, atom.Label)
			}
		}

		p := make([]string, 0)
		products := c.parsedFormulas[c.balancer.separatorPos:]
		for _, prod := range products {
			for _, atom := range prod {
				p = append(p, atom.Label)
			}
		}

		ur := utils.UniqueElems(r)
		up := utils.UniqueElems(p)

		diff := utils.SymmetricDifference(ur, up)

		return diff
	}
	return nil
}

func (c *coeffs) getCoeffs() (MethodResult, error) {
	user := "User"
	nilStr := MethodResult{Method: user, Result: nil}

	diff := c.elementCountValidation()
	if diff != nil {
		return nilStr,
			fmt.Errorf("Cannot balance this reaction, because element(s) %v are only in one part of the reaction", diff)
	}
	coeffs, err := c.calculateCoeffs()
	if err != nil {
		return nilStr, err
	}
	err = c.validateCoeffs(coeffs.Result)
	if err != nil {
		return nilStr, err
	}

	return coeffs, nil
}
