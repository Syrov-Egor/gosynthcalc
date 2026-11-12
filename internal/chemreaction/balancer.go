package chemreaction

import (
	"context"
	"math"
	"sync"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

// balancer computes reaction coefficients. The exact path uses the cached
// rational analysis of the reaction matrix; the legacy path runs the float
// heuristics in [balancingAlgos] and post-processes their results (rounding,
// float validation and optional intification). It is created by
// [ChemicalReaction.Balancer].
type balancer struct {
	reactionMatrix *mat.Dense
	separatorPos   int
	intify         bool
	tolerance      float64
	precision      uint
	bAlgos         *balancingAlgos
	maxDenom       int
	strategy       BalanceStrategy

	exact        *exactMatrix
	analysisOnce sync.Once
	analysisRes  matrixAnalysis
	analysisErr  error
}

// MethodResult is a coefficient list together with the name of the method
// that produced it. Method is [MethodUser] when the coefficients come from the
// reaction string or from [ChemicalReaction.SetCoefficients]. [balancer.Auto]
// reports [MethodInverse] for nullity 1 and [MethodExactLP] for higher
// nullity; the legacy heuristics report [MethodGPInv] and [MethodPPInv].
type MethodResult struct {
	// Method is the name of the algorithm or source that produced Result.
	Method string
	// Result holds one coefficient per compound of the reaction.
	Result []float64
}

// newBalancer returns a balancer for the given reaction matrix, whose first
// separatorPos columns are reactants and the rest products. Precision is the
// rounding precision of legacy float results; tolerance defaults to 1e-8 when
// omitted. intify selects whether legacy results are intified,
// strategy selects exact or legacy routing, and maxDenom (the maximum
// coefficient and fraction denominator, see [balancer.intifyCoefs]) is fixed
// at 1_000_000.
func newBalancer(matrix *mat.Dense, separatorPos int, intify bool, precision uint, strategy BalanceStrategy, tolerance ...float64) *balancer {
	return newBalancerWithExact(matrix, nil, separatorPos, intify, precision, strategy, tolerance...)
}

// newBalancerWithExact is newBalancer with a prebuilt exact reaction matrix.
// When exact is nil the exact matrix is derived from the float matrix on first
// use, which keeps direct test/legacy constructions working.
func newBalancerWithExact(matrix *mat.Dense, exact *exactMatrix, separatorPos int, intify bool, precision uint, strategy BalanceStrategy, tolerance ...float64) *balancer {
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
		strategy:       strategy,
		exact:          exact,
	}
}

// exactMatrixAndAnalysis returns the signed exact matrix and its cached exact
// analysis. The float reaction matrix is converted once per balancer and the
// RREF, rank, nullity and null-space basis are computed once as well.
func (b *balancer) exactMatrixAndAnalysis() (*exactMatrix, matrixAnalysis, error) {
	b.analysisOnce.Do(func() {
		exact := b.exact
		if exact == nil {
			var err error
			exact, err = newExactMatrix(b.reactionMatrix, b.separatorPos)
			if err != nil {
				b.analysisErr = err
				return
			}
			b.exact = exact
		}
		analysis, err := analyze(exact)
		if err != nil {
			b.analysisErr = err
			return
		}
		b.analysisRes = analysis
	})
	return b.exact, b.analysisRes, b.analysisErr
}

// analysis returns the cached exact matrix analysis of the reaction.
func (b *balancer) analysis() (matrixAnalysis, error) {
	_, analysis, err := b.exactMatrixAndAnalysis()
	return analysis, err
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
// result is validated again by [balancer.calculateLegacyByMethod] before it
// is accepted.
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
// within atol. Callers should pass [scaledTolerance] so that large legitimate
// coefficients are not rejected by an absolute threshold.
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

// scaledTolerance returns tolerance scaled by the coefficient magnitude. This
// avoids false failures for large but legitimate coefficients such as
// 1000001 while preserving the configured tolerance for small coefficients.
func scaledTolerance(coefs []float64, tolerance float64) float64 {
	scale := 1.0
	for _, coef := range coefs {
		scale += math.Abs(coef)
	}
	return tolerance * scale
}

// validFloatResult reports whether coefs is an acceptable legacy answer: one
// finite positive coefficient per compound and a weighted reactant/product
// balance within a coefficient-scaled tolerance.
func (b *balancer) validFloatResult(coefs []float64) bool {
	_, matrLength := b.reactionMatrix.Dims()
	if len(coefs) != matrLength || !allPositive(coefs) {
		return false
	}
	for _, coef := range coefs {
		if math.IsNaN(coef) || math.IsInf(coef, 0) {
			return false
		}
	}
	return isReactionBalanced(
		b.bAlgos.ReactantMatrix,
		b.bAlgos.ProductMatrix,
		coefs,
		scaledTolerance(coefs, b.tolerance),
	)
}

// Inv computes the coefficients from the exact null space of the signed
// reaction matrix. It is the exact RREF equivalent of the corrected Thorne
// null-space construction. It is sound and complete when the nullity is 1: the
// unique null direction is checked for zero and mixed-sign entries, normalised
// to a positive primitive integer vector and verified exactly.
//
// Nullity 0, or a null direction with a zero or mixed-sign entry, returns an
// error wrapping [ErrNoPositiveSolution]. Nullity 2 or more returns an error
// wrapping [ErrNullityNotOne]: for higher nullity a null-space basis vector is
// not a coefficient vector, and a positive linear combination must be found
// by the universal solver.
func (b *balancer) Inv() ([]float64, error) {
	return b.calculateExactByMethod(context.Background(), MethodInverse)
}

// GPinv computes the coefficients with the general pseudoinverse algorithm
// of [balancingAlgos.gPInvAlgorithm] (Risteski's method). It is a legacy
// heuristic: its result is validated as a float and may be much larger than a
// primitive integer solution.
func (b *balancer) GPinv() ([]float64, error) {
	return b.calculateLegacyByMethod(context.Background(), MethodGPInv)
}

// PPinv computes the coefficients with the partial pseudoinverse algorithm
// of [balancingAlgos.pPInvAlgorithm] (Risteski's method). It is a legacy
// heuristic: its result is validated as a float and may be much larger than a
// primitive integer solution.
func (b *balancer) PPinv() ([]float64, error) {
	return b.calculateLegacyByMethod(context.Background(), MethodPPInv)
}

// Comb computes the coefficients by brute force with
// [balancingAlgos.combinatorial], enumerating every combination of
// coefficients between 1 and maxCoef. The parallel search is deterministic: it
// returns the lexicographically first valid vector within the bound. ctx
// cancels a running search; an error is returned when no combination balances
// the reaction.
//
// Comb is a bounded, exponential diagnostic. It is incomplete when the true
// solution exceeds maxCoef, and it is not attempted by [balancer.Auto]. Prefer
// [balancer.Auto] or [balancer.UniversalExact] for production balancing; Comb
// remains available for compatibility and comparison.
func (b *balancer) Comb(ctx context.Context, maxCoef uint) ([]float64, error) {
	return b.calculateLegacyByMethod(ctx, MethodComb, maxCoef)
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

// allPositive reports whether every coefficient is strictly greater than 0.
func allPositive(coefs []float64) bool {
	for _, coef := range coefs {
		if coef <= 0 {
			return false
		}
	}
	return true
}
