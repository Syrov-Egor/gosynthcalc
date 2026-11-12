package chemreaction

import (
	"context"
	"fmt"
	"math/big"
)

// defaultSimplexIterationLimit guards against pathological inputs and
// implementation errors. Bland's rule already guarantees finite termination
// for a well-formed Phase-I problem.
const defaultSimplexIterationLimit = 1_000_000

// simplexOptions configures [phaseI]. A zero value uses
// [context.Background] and [defaultSimplexIterationLimit].
type simplexOptions struct {
	ctx           context.Context
	maxIterations int
}

// phaseI solves A c = b with c >= 0 exactly over the rationals. It returns a
// basic feasible solution and true when the system is feasible, nil and false
// when it is infeasible, and an error for malformed input, cancellation or an
// iteration-limit violation.
//
// The Phase-I linear program is
//
//	minimize     sum(a_i)
//	subject to   A c + a = b,
//	             c >= 0, a >= 0.
//
// Rows with a negative right-hand side are negated first, so the artificial
// variables form an initial basis. Bland's rule (smallest negative reduced
// cost for entering, smallest basis index for ratio ties) prevents cycling.
func phaseI(A [][]*big.Rat, b []*big.Rat, options ...simplexOptions) ([]*big.Rat, bool, error) {
	opts := simplexOptions{}
	if len(options) > 0 {
		opts = options[0]
	}
	if opts.ctx == nil {
		opts.ctx = context.Background()
	}
	if opts.maxIterations <= 0 {
		opts.maxIterations = defaultSimplexIterationLimit
	}

	m := len(A)
	if m == 0 {
		return []*big.Rat{}, true, nil
	}
	if len(b) != m {
		return nil, false, fmt.Errorf("%w: phase I has %d rows but %d right-hand sides",
			ErrInvalidMatrix, m, len(b))
	}
	n := len(A[0])
	for i, row := range A {
		if len(row) != n {
			return nil, false, fmt.Errorf("%w: phase I row %d has %d columns, want %d",
				ErrInvalidMatrix, i, len(row), n)
		}
		for j, v := range row {
			if v == nil {
				return nil, false, fmt.Errorf("%w: phase I entry (%d,%d) is nil",
					ErrInvalidMatrix, i, j)
			}
		}
	}
	for i, v := range b {
		if v == nil {
			return nil, false, fmt.Errorf("%w: phase I right-hand side %d is nil",
				ErrInvalidMatrix, i)
		}
	}

	tableau := make([][]*big.Rat, m)
	rhs := make([]*big.Rat, m)
	basis := make([]int, m)

	for i := 0; i < m; i++ {
		tableau[i] = make([]*big.Rat, n+m)
		for j := 0; j < n+m; j++ {
			tableau[i][j] = new(big.Rat)
		}
		for j := 0; j < n; j++ {
			tableau[i][j].Set(A[i][j])
		}

		rhs[i] = new(big.Rat).Set(b[i])
		if rhs[i].Sign() < 0 {
			for j := 0; j < n; j++ {
				tableau[i][j].Neg(tableau[i][j])
			}
			rhs[i].Neg(rhs[i])
		}

		tableau[i][n+i].SetInt64(1)
		basis[i] = n + i
	}

	// Objective w = sum artificials = w0 + sum_j r[j] x_j.
	reduced := make([]*big.Rat, n+m)
	for j := 0; j < n+m; j++ {
		reduced[j] = new(big.Rat)
	}
	w0 := new(big.Rat)
	for i := 0; i < m; i++ {
		w0.Add(w0, rhs[i])
		for j := 0; j < n; j++ {
			reduced[j].Sub(reduced[j], tableau[i][j])
		}
	}

	pivot := func(p, q int) {
		pv := new(big.Rat).Set(tableau[p][q])

		for j := 0; j < n+m; j++ {
			tableau[p][j].Quo(tableau[p][j], pv)
		}
		rhs[p].Quo(rhs[p], pv)

		for i := 0; i < m; i++ {
			if i == p || tableau[i][q].Sign() == 0 {
				continue
			}
			f := new(big.Rat).Set(tableau[i][q])
			for j := 0; j < n+m; j++ {
				tableau[i][j].Sub(tableau[i][j], new(big.Rat).Mul(f, tableau[p][j]))
			}
			rhs[i].Sub(rhs[i], new(big.Rat).Mul(f, rhs[p]))
		}

		if reduced[q].Sign() != 0 {
			f := new(big.Rat).Set(reduced[q])
			for j := 0; j < n+m; j++ {
				reduced[j].Sub(reduced[j], new(big.Rat).Mul(f, tableau[p][j]))
			}
			w0.Add(w0, new(big.Rat).Mul(f, rhs[p]))
		}
		basis[p] = q
	}

	iterations := 0
	for {
		if err := opts.ctx.Err(); err != nil {
			return nil, false, err
		}

		q := -1
		for j := 0; j < n+m; j++ {
			if reduced[j].Sign() < 0 {
				q = j
				break
			}
		}
		if q < 0 {
			break
		}
		if iterations >= opts.maxIterations {
			return nil, false, ErrSimplexIterationLimit
		}

		p := -1
		var best *big.Rat
		for i := 0; i < m; i++ {
			if tableau[i][q].Sign() <= 0 {
				continue
			}
			ratio := new(big.Rat).Quo(rhs[i], tableau[i][q])
			if p < 0 || ratio.Cmp(best) < 0 ||
				(ratio.Cmp(best) == 0 && basis[i] < basis[p]) {
				p = i
				best = ratio
			}
		}
		if p < 0 {
			return nil, false, fmt.Errorf("%w: phase I is unbounded", ErrInvalidMatrix)
		}

		pivot(p, q)
		iterations++
	}

	if w0.Sign() != 0 {
		return nil, false, nil
	}

	c := make([]*big.Rat, n)
	for j := range c {
		c[j] = new(big.Rat)
	}
	for i := 0; i < m; i++ {
		if basis[i] < n {
			c[basis[i]].Set(rhs[i])
		}
	}
	return c, true, nil
}

// balanceExact returns a strictly positive rational vector x with A x = 0. It
// applies the shift c = x - 1, solves the standard-form feasibility problem
// A c = -A 1 with c >= 0 by [phaseI], and returns x = c + 1.
//
// An infeasible system returns an error wrapping [ErrNoPositiveSolution].
// Cancellation and iteration-limit failures are returned unchanged so callers
// can distinguish them with [errors.Is].
func balanceExact(A *exactMatrix, options ...simplexOptions) ([]*big.Rat, error) {
	if err := A.validate(); err != nil {
		return nil, err
	}
	n := A.cols

	if A.rows == 0 {
		x := make([]*big.Rat, n)
		for j := range x {
			x[j] = big.NewRat(1, 1)
		}
		return x, nil
	}

	b := make([]*big.Rat, A.rows)
	for i := 0; i < A.rows; i++ {
		sum := new(big.Rat)
		for j := 0; j < n; j++ {
			sum.Add(sum, A.data[i][j])
		}
		b[i] = new(big.Rat).Neg(sum)
	}

	c, feasible, err := phaseI(A.data, b, options...)
	if err != nil {
		return nil, err
	}
	if !feasible {
		return nil, fmt.Errorf("%w: phase I is infeasible", ErrNoPositiveSolution)
	}

	x := make([]*big.Rat, n)
	for j := 0; j < n; j++ {
		x[j] = new(big.Rat).Add(c[j], big.NewRat(1, 1))
		if x[j].Sign() <= 0 {
			return nil, fmt.Errorf("%w: shift produced a non-positive coefficient", ErrNoPositiveSolution)
		}
	}

	for i := 0; i < A.rows; i++ {
		sum := new(big.Rat)
		for j := 0; j < n; j++ {
			sum.Add(sum, new(big.Rat).Mul(A.data[i][j], x[j]))
		}
		if sum.Sign() != 0 {
			return nil, fmt.Errorf("%w: phase I produced a non-balancing vector", ErrInvalidMatrix)
		}
	}
	return x, nil
}

// universalExact runs the exact Phase-I simplex and returns primitive positive
// integer coefficients verified against A.
func (b *balancer) universalExact(ctx context.Context) ([]*big.Int, error) {
	A, _, err := b.exactMatrixAndAnalysis()
	if err != nil {
		return nil, err
	}
	x, err := balanceExact(A, simplexOptions{ctx: ctx})
	if err != nil {
		return nil, err
	}
	z, err := primitiveIntegers(x)
	if err != nil {
		return nil, err
	}
	if !validExactResult(A, z) {
		return nil, fmt.Errorf("exact simplex solution failed validation: %w", ErrInvalidMatrix)
	}
	return z, nil
}

// UniversalExact balances the reaction with the exact rational Phase-I
// simplex and returns primitive positive integer coefficients. It is complete
// for every rational reaction matrix; an infeasible reaction returns an error
// wrapping [ErrNoPositiveSolution].
func (b *balancer) UniversalExact(ctx context.Context) ([]*big.Int, error) {
	return b.universalExact(ctx)
}

// Universal is the float64 convenience wrapper around [balancer.UniversalExact].
// It returns [ErrCoefficientOverflow] if a primitive coefficient cannot be
// represented exactly as a float64.
func (b *balancer) Universal(ctx context.Context) ([]float64, error) {
	z, err := b.UniversalExact(ctx)
	if err != nil {
		return nil, err
	}
	return bigIntsToFloat64(z)
}

// solveSimplexExact runs the exact simplex for a reaction whose analysis is
// already known. Nullity 0 is rejected before starting the LP.
func (b *balancer) solveSimplexExact(ctx context.Context, a matrixAnalysis) ([]*big.Int, error) {
	if a.nullity == 0 {
		return nil, fmt.Errorf("%w: reaction matrix has trivial null space", ErrNoPositiveSolution)
	}
	return b.universalExact(ctx)
}

// solveSimplex is the float64 convenience wrapper around
// [balancer.solveSimplexExact].
func (b *balancer) solveSimplex(ctx context.Context, a matrixAnalysis) ([]float64, error) {
	z, err := b.solveSimplexExact(ctx, a)
	if err != nil {
		return nil, err
	}
	return bigIntsToFloat64(z)
}
