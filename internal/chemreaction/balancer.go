package chemreaction

import (
	"context"
	"fmt"
	"math"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

// balancer computes reaction coefficients with the matrix algorithms of
// [balancingAlgos] and post-processes their results: it rounds them, checks
// that they really balance the reaction and, when intify is set, converts
// them to the smallest integer equivalents. It is created by [ChemicalReaction.Balancer].
type balancer struct {
	reactionMatrix *mat.Dense
	separatorPos   int
	intify         bool
	tolerance      float64
	precision      uint
	bAlgos         *balancingAlgos
	maxDenom       int
}

// MethodResult is a coefficient list together with the name of the method
// that produced it. Method is "User" when the coefficients come from the
// reaction string or from [ChemicalReaction.SetCoefficients], and
// "inverse", "general pseudoinverse" or "partial pseudoinverse" when
// they come from [balancer.Auto].
type MethodResult struct {
	// Method is the name of the algorithm or source that produced Result.
	Method string
	// Result holds one coefficient per compound of the reaction.
	Result []float64
}

// newBalancer returns a balancer for the given reaction matrix, whose first
// separatorPos columns are reactants and the rest products. Precision is the
// rounding precision of the coefficients; tolerance defaults to 1e-8 when
// omitted. intify selects whether integer coefficients are preferred, and
// maxDenom (the maximum coefficient and fraction denominator, see
// [balancer.intifyCoefs]) is fixed at 1_000_000.
func newBalancer(matrix *mat.Dense, separatorPos int, intify bool, precision uint, tolerance ...float64) *balancer {
	var tol float64
	if tolerance == nil {
		tol = 1e-8
	} else {
		tol = tolerance[0]
	}

	bAlgos := newBalancingAlgos(matrix, separatorPos, tol)

	return &balancer{
		reactionMatrix: matrix,
		separatorPos:   separatorPos,
		intify:         intify,
		tolerance:      tol,
		precision:      precision,
		bAlgos:         bAlgos,
		maxDenom:       1_000_000,
	}
}

// reactionRankAndNullity returns the rank and the nullity (number of columns
// minus rank) of the reaction matrix. The algorithms in [balancingAlgos]
// receive these values from the balancer instead of calculating them again.
func (b *balancer) reactionRankAndNullity() (rank, nullity int, err error) {
	rank, _, err = matrixRank(b.reactionMatrix, b.tolerance)
	if err != nil {
		return 0, 0, fmt.Errorf("error computing rank of the reaction matrix: %w", err)
	}
	_, cols := b.reactionMatrix.Dims()
	return rank, cols - rank, nil
}

// reactantRank returns the rank of the reactant half of the reaction matrix.
// It is needed by [balancingAlgos.pPInvAlgorithm] together with the rank of
// the whole matrix.
func (b *balancer) reactantRank() (int, error) {
	rank, _, err := matrixRank(b.bAlgos.ReactantMatrix, b.tolerance)
	if err != nil {
		return 0, fmt.Errorf("error computing rank of the reactant matrix: %w", err)
	}
	return rank, nil
}

// intifyCoefs converts fractional coefficients into the smallest whole
// numbers representing the same ratios: every coefficient is turned into a
// fraction whose denominator is bounded by maxDenom, the denominators are
// cleared with their least common multiple, and the resulting integers are
// divided by their greatest common divisor.
//
// The conversion is all-or-nothing: if any coefficient is not finite or
// exceeds limit, cannot be represented as an int64 fraction, would make the
// LCM overflow or exceed 1e15, flips sign along the way or ends up negative
// or oversized, the input coefficients are returned unchanged so that
// callers keep working with correct floats rather than wrong integers. The
// result is validated again by [balancer.calculateByMethod] before it is
// accepted.
func (b *balancer) intifyCoefs(coefs []float64, limit int) []float64 {
	initialCoefficients := make([]float64, len(coefs))
	copy(initialCoefficients, coefs)

	fractions := make([]utils.SimpleFraction, len(coefs))
	denominators := make([]int64, len(coefs))

	for i, coef := range coefs {
		if math.IsNaN(coef) || math.IsInf(coef, 0) || math.Abs(coef) > float64(limit) {
			return initialCoefficients
		}
		frac := utils.NewSimpleFraction(coef, int64(b.maxDenom))
		if frac.Den == 0 {
			return initialCoefficients
		}
		fractions[i] = frac
		denominators[i] = frac.Den
	}

	lcm := utils.FindLCMSliceInt64(denominators)
	if lcm < 0 || lcm > 1e15 {
		return initialCoefficients
	}

	vals := make([]int64, len(fractions))
	for i, frac := range fractions {
		vals[i] = frac.Num * (lcm / frac.Den)

		if vals[i] < 0 && frac.Num > 0 {
			return initialCoefficients
		}
	}

	gcd := utils.FindGCDSliceInt64(vals)
	if gcd == 0 {
		return initialCoefficients
	}

	coefficients := make([]int64, len(vals))
	for i, val := range vals {
		coefficients[i] = val / gcd
	}

	// Chemical coefficients must stay positive, so a negative value here
	// means the input was not a valid reaction: fall back to the validated
	// floats instead of silently flipping the sign, which would forge a
	// plausible-looking answer out of a wrong one. The same fallback catches
	// values that clearing denominators pushed past limit.
	for _, coeff := range coefficients {
		if coeff < 0 || coeff > int64(limit) {
			return initialCoefficients
		}
	}

	result := make([]float64, len(coefficients))
	for i, coeff := range coefficients {
		result[i] = float64(coeff)
	}

	return result
}

// isReactionBalanced reports whether coefficients balance the reaction: the
// reactant matrix and the product matrix are each multiplied by their half
// of the coefficient vector, and the two weighted column sums must be equal
// within atol.
func isReactionBalanced(reactantMatrix *mat.Dense, productMatrix *mat.Dense, coefs []float64, atol float64) bool {
	reactantRows, reactantCols := reactantMatrix.Dims()
	productRows, productCols := productMatrix.Dims()
	separatorPos := reactantCols
	reactantCoefs := coefs[:separatorPos]
	productCoefs := coefs[separatorPos:]
	reacSum := make([]float64, reactantRows)
	prodSum := make([]float64, productRows)
	mulAndSumFl(reactantMatrix, reactantCoefs, reacSum, reactantRows, reactantCols)
	mulAndSumFl(productMatrix, productCoefs, prodSum, productRows, productCols)

	return floats.EqualApprox(reacSum, prodSum, atol)
}

// calculateByMethod runs one of the named algorithms — "inv", "gpinv",
// "ppinv" or "comb" (only "comb" reads maxCoef, the largest coefficient the
// search may try) — and post-processes its result: the coefficients are
// rounded to precision+2 decimals and must pass [balancer.validResult],
// meaning their count, positivity and balance are all verified within
// tolerance. If intify is enabled, [balancer.intifyCoefs] is attempted on top
// and is only kept when it still passes the same validation.
//
// Any failure is reported as an error ("can't balance reaction by %s
// method", "wrong coefficients" or "no method %s"). ctx is used to cancel
// a running "comb" search.
func (b *balancer) calculateByMethod(ctx context.Context, method string, maxCoef ...uint) ([]float64, error) {
	var coefficients []float64
	var err error
	errm := fmt.Errorf("can't balance reaction by %s method", method)

	switch method {
	case "inv":
		rank, nullity, rankErr := b.reactionRankAndNullity()
		if rankErr != nil {
			return nil, errm
		}
		coefficients, err = b.bAlgos.invAlgorithm(rank, nullity)
		if err != nil {
			return nil, errm
		}
	case "gpinv":
		rank, nullity, rankErr := b.reactionRankAndNullity()
		if rankErr != nil {
			return nil, errm
		}
		coefficients, err = b.bAlgos.gPInvAlgorithm(rank, nullity)
		if err != nil {
			return nil, errm
		}
	case "ppinv":
		rank, nullity, rankErr := b.reactionRankAndNullity()
		if rankErr != nil {
			return nil, errm
		}
		reactantRank, rankErr := b.reactantRank()
		if rankErr != nil {
			return nil, errm
		}
		coefficients, err = b.bAlgos.pPInvAlgorithm(rank, nullity, reactantRank)
		if err != nil {
			return nil, errm
		}
	case "comb":
		coefficients = b.bAlgos.combinatorial(ctx, maxCoef[0])
		if coefficients == nil {
			return nil, errm
		}
	default:
		return nil, fmt.Errorf("no method %s", method)
	}

	return b.finalizeCoefficients(coefficients)
}

// finalizeCoefficients rounds a raw algorithm result, validates it and, when
// intify is enabled, replaces it with the smallest integer coefficients that
// still pass the same validation.
func (b *balancer) finalizeCoefficients(coefficients []float64) ([]float64, error) {
	coefficients = utils.RoundFloatS(coefficients, b.precision+2)
	_, matrLength := b.reactionMatrix.Dims()

	if !b.validResult(coefficients, matrLength) {
		return nil, fmt.Errorf("wrong coefficients")
	}

	if b.intify {
		intified := b.intifyCoefs(coefficients, b.maxDenom)
		if b.validResult(intified, matrLength) {
			coefficients = intified
		}
	}
	return coefficients, nil
}

// validResult reports whether coefs is an acceptable answer: it must contain
// one coefficient per compound, every coefficient must be positive, and the
// reaction must be balanced within tolerance.
func (b *balancer) validResult(coefs []float64, matrLength int) bool {
	return len(coefs) == matrLength &&
		utils.AllPositive(coefs) &&
		isReactionBalanced(
			b.bAlgos.ReactantMatrix,
			b.bAlgos.ProductMatrix,
			coefs,
			b.tolerance,
		)
}

// Inv computes the coefficients with the matrix inverse algorithm of
// [balancingAlgos.invAlgorithm] (Thorne's method).
func (b *balancer) Inv() ([]float64, error) {
	res, err := b.calculateByMethod(context.Background(), "inv")
	if err != nil {
		return nil, err
	}
	return res, nil
}

// GPinv computes the coefficients with the general pseudoinverse algorithm
// of [balancingAlgos.gPInvAlgorithm] (Risteski's method).
func (b *balancer) GPinv() ([]float64, error) {
	res, err := b.calculateByMethod(context.Background(), "gpinv")
	if err != nil {
		return nil, err
	}
	return res, nil
}

// PPinv computes the coefficients with the partial pseudoinverse algorithm
// of [balancingAlgos.pPInvAlgorithm] (Risteski's method).
func (b *balancer) PPinv() ([]float64, error) {
	res, err := b.calculateByMethod(context.Background(), "ppinv")
	if err != nil {
		return nil, err
	}
	return res, nil
}

// Comb computes the coefficients by brute force with
// [balancingAlgos.combinatorial], enumerating every combination of
// coefficients between 1 and maxCoef. ctx cancels a running search; an error
// is returned when no combination balances the reaction.
//
// Unlike the matrix methods, this one is not attempted by [balancer.Auto].
func (b *balancer) Comb(ctx context.Context, maxCoef uint) ([]float64, error) {
	res, err := b.calculateByMethod(ctx, "comb", maxCoef)
	if err != nil {
		return nil, err
	}
	return res, nil
}

// Auto chooses a balancing algorithm from the rank and nullity of the
// reaction matrix:
//
//   - nullity 0: no balanced equation exists, so an error is returned;
//   - nullity 1: the unique (up to scale) solution is computed with
//     [balancingAlgos.invAlgorithm] ("inverse");
//   - nullity > 1: [balancingAlgos.gPInvAlgorithm] ("general pseudoinverse")
//     and [balancingAlgos.pPInvAlgorithm] ("partial pseudoinverse") are both
//     tried, and the valid result with the lowest sum of raw coefficients is
//     returned. The sums are compared before intification; the chosen
//     coefficients then go through the usual rounding, validation and
//     intification.
//
// If no algorithm produces valid coefficients it reports "can't balance this
// reaction by any method".
func (b *balancer) Auto() (MethodResult, error) {
	rank, nullity, err := b.reactionRankAndNullity()
	if err != nil {
		return MethodResult{Method: "", Result: nil},
			fmt.Errorf("can't balance this reaction by any method: %w", err)
	}

	if nullity == 0 {
		return MethodResult{Method: "", Result: nil},
			fmt.Errorf("the reaction cannot be balanced (nullity 0)")
	}

	if nullity == 1 {
		coefficients, err := b.bAlgos.invAlgorithm(rank, nullity)
		if err != nil {
			return MethodResult{Method: "", Result: nil},
				fmt.Errorf("can't balance this reaction by any method: %w", err)
		}
		coefficients, err = b.finalizeCoefficients(coefficients)
		if err != nil {
			return MethodResult{Method: "", Result: nil},
				fmt.Errorf("can't balance this reaction by any method: %w", err)
		}
		return MethodResult{Method: "inverse", Result: coefficients}, nil
	}

	reactantRank, err := b.reactantRank()
	if err != nil {
		return MethodResult{Method: "", Result: nil},
			fmt.Errorf("can't balance this reaction by any method: %w", err)
	}

	type candidate struct {
		method string
		coefs  []float64
		sum    float64
	}
	candidates := make([]candidate, 0, 2)

	if raw, err := b.bAlgos.gPInvAlgorithm(rank, nullity); err == nil {
		if coefs, err := b.finalizeCoefficients(raw); err == nil {
			candidates = append(candidates, candidate{
				method: "general pseudoinverse",
				coefs:  coefs,
				sum:    utils.SumAbs(raw),
			})
		}
	}
	if raw, err := b.bAlgos.pPInvAlgorithm(rank, nullity, reactantRank); err == nil {
		if coefs, err := b.finalizeCoefficients(raw); err == nil {
			candidates = append(candidates, candidate{
				method: "partial pseudoinverse",
				coefs:  coefs,
				sum:    utils.SumAbs(raw),
			})
		}
	}

	if len(candidates) == 0 {
		return MethodResult{Method: "", Result: nil},
			fmt.Errorf("can't balance this reaction by any method")
	}

	best := candidates[0]
	for _, cand := range candidates[1:] {
		if cand.sum < best.sum {
			best = cand
		}
	}
	return MethodResult{Method: best.method, Result: best.coefs}, nil
}

// mulAndSumFl computes result = matrix · vector in place for a floating
// point coefficient vector (rows weighted sums of the matrix columns).
func mulAndSumFl(matrix *mat.Dense, vector []float64, result []float64, rows int, cols int) {
	for row := range rows {
		result[row] = 0
		for col := range cols {
			result[row] += matrix.At(row, col) * vector[col]
		}
	}
}
