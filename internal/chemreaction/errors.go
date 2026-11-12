package chemreaction

import "errors"

// Sentinel errors used by the exact balancing pipeline. They are wrapped with
// context and can be tested with [errors.Is].
var (
	// ErrNoPositiveSolution means the reaction has no strictly positive
	// coefficient vector.
	ErrNoPositiveSolution = errors.New("reaction has no strictly positive solution")
	// ErrNoNullspace means the reaction matrix has a trivial null space.
	ErrNoNullspace = errors.New("reaction matrix has trivial null space")
	// ErrNullityNotOne means the requested operation only applies to a
	// matrix whose nullity is exactly one.
	ErrNullityNotOne = errors.New("nullity is not one")
	// ErrMethodNotPositive means a heuristic method produced a candidate
	// that is not strictly positive. It does not prove infeasibility.
	ErrMethodNotPositive = errors.New("method produced a non-positive candidate")
	// ErrSimplexIterationLimit means the exact simplex exceeded its safety
	// iteration limit.
	ErrSimplexIterationLimit = errors.New("simplex iteration limit exceeded")
	// ErrCoefficientOverflow means a coefficient cannot be represented
	// exactly as a float64.
	ErrCoefficientOverflow = errors.New("coefficient exceeds exact float64 range")
	// ErrNonFiniteMatrix means the reaction matrix contains NaN or an
	// infinity.
	ErrNonFiniteMatrix = errors.New("reaction matrix contains non-finite values")
	// ErrInvalidMatrix means the matrix, its dimensions or a helper input is
	// invalid.
	ErrInvalidMatrix = errors.New("invalid reaction matrix")
)
