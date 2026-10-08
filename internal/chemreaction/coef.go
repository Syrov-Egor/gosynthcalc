package chemreaction

import (
	"fmt"

	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

// coeffs computes and validates the coefficients of a reaction according to
// the selected [Mode]. Element-set validation, coefficient calculation for
// the mode, coefficient list validation.
type coeffs struct {
	mode               Mode
	parsedFormulas     [][]chemformula.Atom
	decomposedReaction *reactionDecomposer
	balancer           *balancer
}

// calculateCoeffs produces the coefficients for the configured mode:
//
//   - [Force] takes the coefficients written in the reaction string as they
//     are, without any balance check (Method "User");
//   - [Check] takes them only if they balance the reaction, and fails with
//     "reaction is not balanced" otherwise (Method "User");
//   - [Balance] computes them automatically with [balancer.Auto] and reports
//     the algorithm that succeeded.
//
// Any other mode fails with "no such mode".
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
			return MethodResult{Method: user, Result: nil}, fmt.Errorf("reaction is not balanced")
		}

	case Balance:
		coefs, err := c.balancer.Auto()
		if err != nil {
			return MethodResult{Method: user, Result: nil}, err
		}
		return coefs, nil

	default:
		return MethodResult{Method: user, Result: nil}, fmt.Errorf("no such mode %d", c.mode)
	}
}

// validateCoeffs checks that every coefficient is strictly positive and that
// there is exactly one per compound of the reaction.
func (c *coeffs) validateCoeffs(coefs []float64) error {
	_, cols := c.balancer.reactionMatrix.Dims()

	switch {
	case !allPositive(coefs):
		return fmt.Errorf("some coefs in %v are negative or 0", coefs)
	case len(coefs) != cols:
		return fmt.Errorf("number of coefs should be equal %d, got %d", cols, len(coefs))
	default:
		return nil
	}
}

// elementCountValidation returns the [utils.SymmetricDifference] of the
// element sets of the reactant and product parts: a non-empty result means
// some element occurs on only one side of the reaction, which no
// coefficients can fix, and the caller turns it into an error. In [Force]
// mode the check is skipped and nil is returned, since the user takes
// responsibility for the reaction.
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

// getCoeffs runs the whole pipeline — element-set validation, mode-specific
// calculation and list validation — and returns the coefficients together
// with the name of the method that produced them.
func (c *coeffs) getCoeffs() (MethodResult, error) {
	user := "User"
	nilStr := MethodResult{Method: user, Result: nil}

	diff := c.elementCountValidation()
	if diff != nil {
		return nilStr,
			fmt.Errorf("cannot balance this reaction, because element(s) %v are only in one part of the reaction", diff)
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
