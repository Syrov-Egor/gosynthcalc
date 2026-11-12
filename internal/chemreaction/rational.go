package chemreaction

import (
	"fmt"
	"math"
	"math/big"
	"strconv"
	"strings"
)

// ratFromFloat converts a float64 to an exact rational number.
//
// The conversion goes through the shortest fixed-point decimal that parses
// back to v. Formula atom counts no longer rely on this path: the formula
// parser carries exact rational amounts in parallel and the balancer uses
// them directly. This helper remains for other float inputs, such as
// user-supplied coefficients or legacy results.
//
// NaN and infinities are rejected with [ErrNonFiniteMatrix].
func ratFromFloat(v float64) (*big.Rat, error) {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return nil, fmt.Errorf("%w: %v", ErrNonFiniteMatrix, v)
	}
	s := strconv.FormatFloat(v, 'f', -1, 64)
	r, ok := new(big.Rat).SetString(s)
	if !ok {
		return nil, fmt.Errorf("%w: cannot convert %q to a rational", ErrInvalidMatrix, s)
	}
	return r, nil
}

// ratFromDecimal parses an exact rational number from its decimal string
// form. Surrounding whitespace is ignored. Malformed input and the special
// values "NaN" and "Inf" are rejected with [ErrInvalidMatrix].
func ratFromDecimal(s string) (*big.Rat, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return nil, fmt.Errorf("%w: empty rational string", ErrInvalidMatrix)
	}
	r, ok := new(big.Rat).SetString(s)
	if !ok {
		return nil, fmt.Errorf("%w: invalid rational %q", ErrInvalidMatrix, s)
	}
	return r, nil
}

// bigIntsToFloat64 converts a primitive integer vector to float64, checking
// that every value is represented exactly. If a value would lose precision or
// overflow, [ErrCoefficientOverflow] is returned instead of a silently
// corrupted coefficient.
func bigIntsToFloat64(z []*big.Int) ([]float64, error) {
	out := make([]float64, len(z))
	for i, v := range z {
		if v == nil {
			return nil, fmt.Errorf("%w: nil coefficient at position %d", ErrInvalidMatrix, i)
		}
		f, _ := new(big.Float).SetInt(v).Float64()
		if math.IsInf(f, 0) || math.IsNaN(f) {
			return nil, fmt.Errorf("%w: %s", ErrCoefficientOverflow, v.String())
		}
		back := new(big.Rat).SetInt(v)
		roundTrip := new(big.Rat).SetFloat64(f)
		if roundTrip == nil || roundTrip.Cmp(back) != 0 {
			return nil, fmt.Errorf("%w: %s", ErrCoefficientOverflow, v.String())
		}
		out[i] = f
	}
	return out, nil
}

// exactBalances reports whether A z = 0 exactly, with z interpreted as a
// rational vector. It returns false for a malformed matrix or a vector of the
// wrong length.
func exactBalances(A *exactMatrix, z []*big.Int) bool {
	if A == nil || A.cols <= 0 || A.rows < 0 || len(z) != A.cols {
		return false
	}
	if len(A.data) != A.rows {
		return false
	}
	for i := 0; i < A.rows; i++ {
		if len(A.data[i]) != A.cols {
			return false
		}
		sum := new(big.Rat)
		for j := 0; j < A.cols; j++ {
			if z[j] == nil {
				return false
			}
			term := new(big.Rat).Mul(A.data[i][j], new(big.Rat).SetInt(z[j]))
			sum.Add(sum, term)
		}
		if sum.Sign() != 0 {
			return false
		}
	}
	return true
}

// exactIntegerCoefficients converts a float coefficient vector to exact
// integers when every value is finite and integral. It returns false when any
// value is fractional or non-finite. Values are interpreted through their
// shortest decimal representation, consistently with [ratFromFloat].
func exactIntegerCoefficients(coefs []float64) ([]*big.Int, bool) {
	if len(coefs) == 0 {
		return nil, false
	}
	out := make([]*big.Int, len(coefs))
	for i, f := range coefs {
		if math.IsNaN(f) || math.IsInf(f, 0) {
			return nil, false
		}
		r, err := ratFromFloat(f)
		if err != nil || !r.IsInt() {
			return nil, false
		}
		out[i] = new(big.Int).Set(r.Num())
	}
	return out, true
}

// validExactResult reports whether z is an acceptable exact answer for A: one
// strictly positive integer per compound and A z = 0 exactly.
func validExactResult(A *exactMatrix, z []*big.Int) bool {
	if A == nil || len(z) != A.cols {
		return false
	}
	for _, v := range z {
		if v == nil || v.Sign() <= 0 {
			return false
		}
	}
	return exactBalances(A, z)
}

// primitiveIntegers clears the denominators of a strictly positive rational
// vector and divides the resulting integers by their greatest common divisor.
// The result is the unique primitive positive integer vector on the same ray.
// Empty, non-positive or malformed input is rejected.
func primitiveIntegers(x []*big.Rat) ([]*big.Int, error) {
	if len(x) == 0 {
		return nil, fmt.Errorf("%w: empty coefficient vector", ErrInvalidMatrix)
	}

	lcm := big.NewInt(1)
	for i, v := range x {
		if v == nil {
			return nil, fmt.Errorf("%w: nil coefficient at position %d", ErrInvalidMatrix, i)
		}
		if v.Sign() <= 0 {
			return nil, fmt.Errorf("%w: coefficient at position %d is not positive", ErrNoPositiveSolution, i)
		}
		den := new(big.Int).Set(v.Denom())
		g := new(big.Int).GCD(nil, nil, lcm, den)
		lcm.Quo(lcm, g)
		lcm.Mul(lcm, den)
	}

	nums := make([]*big.Int, len(x))
	gcd := new(big.Int)
	for i, v := range x {
		factor := new(big.Int).Quo(new(big.Int).Set(lcm), v.Denom())
		nums[i] = new(big.Int).Mul(v.Num(), factor)
		gcd.GCD(nil, nil, gcd, new(big.Int).Abs(nums[i]))
	}
	if gcd.Sign() == 0 {
		return nil, fmt.Errorf("%w: zero coefficient vector", ErrNoPositiveSolution)
	}
	for i := range nums {
		nums[i].Quo(nums[i], gcd)
	}
	return nums, nil
}
