package chemreaction

import (
	"fmt"
	"math/big"
)

// nullityOneSolution extracts the unique null direction of a nullity-one
// analysis and checks whether it can be scaled to a strictly positive
// coefficient vector.
//
// A nullity-one null space is a line span{v}. A strictly positive solution
// exists exactly when v has no zero entries and all its entries have the same
// sign. If v is strictly negative it is negated; if it has a zero or mixed
// signs there is no strictly positive solution.
//
// Nullity 0 is reported as [ErrNoPositiveSolution]; nullity 2 or more as
// [ErrNullityNotOne]. The returned rational vector is not yet cleared to
// primitive integers.
func nullityOneSolution(a matrixAnalysis) ([]*big.Rat, error) {
	switch {
	case a.nullity == 0:
		return nil, fmt.Errorf("%w: reaction matrix has trivial null space", ErrNoPositiveSolution)
	case a.nullity != 1:
		return nil, fmt.Errorf("%w: got nullity %d", ErrNullityNotOne, a.nullity)
	}
	if len(a.basis) == 0 || len(a.basis[0]) != 1 {
		return nil, fmt.Errorf("%w: malformed null-space basis", ErrInvalidMatrix)
	}

	v := make([]*big.Rat, len(a.basis))
	negatives := 0
	zeros := 0
	for i, row := range a.basis {
		if len(row) != 1 || row[0] == nil {
			return nil, fmt.Errorf("%w: malformed null-space basis row %d", ErrInvalidMatrix, i)
		}
		v[i] = new(big.Rat).Set(row[0])
		switch v[i].Sign() {
		case -1:
			negatives++
		case 0:
			zeros++
		}
	}

	if zeros > 0 {
		return nil, fmt.Errorf("%w: null direction has %d zero entries",
			ErrNoPositiveSolution, zeros)
	}
	if negatives == len(v) {
		for i := range v {
			v[i].Neg(v[i])
		}
		return v, nil
	}
	if negatives > 0 {
		return nil, fmt.Errorf("%w: null direction has mixed signs", ErrNoPositiveSolution)
	}
	return v, nil
}

// solveNullityOneExact returns the primitive positive integer coefficients of
// a nullity-one reaction, verified against the exact signed matrix. It is the
// exact core of [balancer.Inv].
func (b *balancer) solveNullityOneExact(a matrixAnalysis) ([]*big.Int, error) {
	A, _, err := b.exactMatrixAndAnalysis()
	if err != nil {
		return nil, err
	}

	x, err := nullityOneSolution(a)
	if err != nil {
		return nil, err
	}
	z, err := primitiveIntegers(x)
	if err != nil {
		return nil, err
	}
	if !validExactResult(A, z) {
		return nil, fmt.Errorf("nullity-one solution failed exact validation: %w", ErrInvalidMatrix)
	}
	return z, nil
}

// solveNullityOne returns the primitive integer coefficients of a nullity-one
// reaction as float64 values. It returns [ErrCoefficientOverflow] when a
// coefficient cannot be represented exactly as a float64.
func (b *balancer) solveNullityOne(a matrixAnalysis) ([]float64, error) {
	z, err := b.solveNullityOneExact(a)
	if err != nil {
		return nil, err
	}
	return bigIntsToFloat64(z)
}
