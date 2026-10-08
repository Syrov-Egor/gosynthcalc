// Package utils provides small generic helpers shared by the other
// gosynthcalc packages: counting and deduplicating string elements, summing
// and rounding floats, approximating floats by fractions and the integer
// LCM/GCD arithmetic used to turn fractional reaction coefficients into
// whole numbers.
package utils

import (
	"math"
	"strings"
)

// StringCounter returns how many times every rune of s occurs, mapping each
// distinct character to its count.
func StringCounter(s string) map[string]int {
	counts := make(map[string]int)
	for _, char := range s {
		counts[string(char)]++
	}
	return counts
}

// UniqueElems returns the distinct elements of atomsList in the order of
// their first occurrence. It is used to build the element sets of the two
// sides of a reaction before comparing them.
func UniqueElems(atomsList []string) []string {
	seen := make(map[string]bool)
	uniqueAtomsList := []string{}

	for _, atom := range atomsList {
		if !seen[atom] {
			seen[atom] = true
			uniqueAtomsList = append(uniqueAtomsList, atom)
		}
	}
	return uniqueAtomsList
}

// SumFloatS returns the sum of all values in s (0 for an empty slice).
func SumFloatS(s []float64) float64 {
	var sum = 0.0
	for _, el := range s {
		sum += el
	}
	return sum
}

// pow10 caches the powers of ten for precisions 0-22 so that rounding does
// not have to call math.Pow in the common cases.
var pow10 = [23]float64{
	1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
	1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22,
}

// roundByRatio rounds val to a single unit of ratio (the number of decimal
// places that ratio stands for), with ties going away from zero.
func roundByRatio(val, ratio float64) float64 {
	return math.Round(val*ratio) / ratio
}

// RoundFloat rounds val to precision decimal places (ties away from zero).
// Precisions up to 22 reuse the precomputed pow10 table; larger precisions
// fall back to math.Pow.
func RoundFloat(val float64, precision uint) float64 {
	if precision < uint(len(pow10)) {
		return roundByRatio(val, pow10[precision])
	}
	return roundByRatio(val, math.Pow(10, float64(precision)))
}

// RoundFloatS rounds every value of s to precision decimal places and returns
// the results as a new slice.
func RoundFloatS(s []float64, precision uint) []float64 {
	res := make([]float64, len(s))
	if precision < uint(len(pow10)) {
		ratio := pow10[precision]
		for i, val := range s {
			res[i] = roundByRatio(val, ratio)
		}
		return res
	}
	for i, val := range s {
		res[i] = RoundFloat(val, precision)
	}
	return res
}

// SimpleFraction is a rational approximation of a float64 value: the
// represented number is Num/Den. A zero Den marks a value that could not be
// represented as an int64 fraction at all, telling callers to fall back to
// their original floating point input instead of trusting a wrong integer.
type SimpleFraction struct {
	Num, Den int64
}

// unrepresentableFraction marks a value that cannot be turned into an int64
// fraction: the zero denominator tells callers to fall back to their original
// (floating point) value instead of trusting a wrong integer.
func unrepresentableFraction() SimpleFraction {
	return SimpleFraction{0, 0}
}

// NewSimpleFraction converts f into a fraction whose denominator never exceeds
// maxDenominator. The value itself is never clamped: exact integers are kept
// as-is and values that cannot be bounded by maxDenominator are returned as
// SimpleFraction{0, 0} or as a nearest-integer approximation, which callers
// have to validate against their own tolerance.
//
// Zero, NaN and the infinities are returned as SimpleFraction{0, 1}.
// Otherwise the sign is stored in Num and the magnitude is searched with a
// bounded Stern-Brocot mediant scan.
func NewSimpleFraction(f float64, maxDenominator int64) SimpleFraction {
	if math.IsInf(f, 0) || math.IsNaN(f) || f == 0 {
		return SimpleFraction{0, 1}
	}

	sign := int64(1)
	if f < 0 {
		sign = -1
		f = -f
	}

	// Exact integers are always representable by denominator 1, no matter how
	// large they are. Clamping them to maxDenominator would silently change the
	// value, so the value itself is preserved here instead.
	if f == math.Floor(f) {
		if f < float64(math.MaxInt64) {
			return SimpleFraction{sign * int64(f), 1}
		}
		return unrepresentableFraction()
	}

	// |f| >= maxDenominator with a fractional part: no fraction bounded by
	// maxDenominator can hold it, so return the nearest integer instead of
	// clamping the value down to maxDenominator. Callers must verify that the
	// approximation still satisfies their constraints.
	if f >= float64(maxDenominator) {
		if f < float64(math.MaxInt64) {
			return SimpleFraction{sign * int64(math.Round(f)), 1}
		}
		return unrepresentableFraction()
	}

	p0, q0 := int64(0), int64(1)
	p1, q1 := int64(1), int64(0)

	for {
		medDen := q0 + q1
		if medDen > maxDenominator {
			break
		}
		medNum := p0 + p1
		medVal := float64(medNum) / float64(medDen)
		switch {
		case medVal < f:
			lo, hi := int64(1), maxDenominator
			for lo < hi {
				mid := lo + (hi-lo+1)/2
				nextDen := q0 + mid*q1
				if nextDen > maxDenominator {
					hi = mid - 1
					continue
				}
				nextNum := p0 + mid*p1
				if float64(nextNum)/float64(nextDen) < f {
					lo = mid
				} else {
					hi = mid - 1
				}
			}
			p0, q0 = p0+lo*p1, q0+lo*q1
		case medVal > f:
			lo, hi := int64(1), maxDenominator
			for lo < hi {
				mid := lo + (hi-lo+1)/2
				nextDen := mid*q0 + q1
				if nextDen > maxDenominator {
					hi = mid - 1
					continue
				}
				nextNum := mid*p0 + p1
				if float64(nextNum)/float64(nextDen) > f {
					lo = mid
				} else {
					hi = mid - 1
				}
			}
			p1, q1 = lo*p0+p1, lo*q0+q1
		default:
			return SimpleFraction{sign * medNum, medDen}
		}
	}
	bestNum, bestDen := p0, q0
	bestError := math.Abs(f - float64(p0)/float64(q0))

	if q1 > 0 && q1 <= maxDenominator {
		upperError := math.Abs(f - float64(p1)/float64(q1))
		if upperError < bestError {
			bestNum, bestDen = p1, q1
		}
	}

	return SimpleFraction{sign * bestNum, bestDen}
}

// FindLCMSliceInt64 returns the least common multiple of nums as a
// non-negative int64; the LCM of an empty slice is 1. If the result would
// overflow int64 or exceed about 1e15 in magnitude (the budget the
// coefficient intifier works with), -1 is returned so callers can bail out
// instead of trusting a wrapped value.
func FindLCMSliceInt64(nums []int64) int64 {
	if len(nums) == 0 {
		return 1
	}

	result := nums[0]
	for i := 1; i < len(nums); i++ {
		if willLCMOverflow(result, nums[i]) {
			return -1
		}
		result = lcmInt64(result, nums[i])
	}
	return result
}

// willLCMOverflow reports whether lcmInt64(a, b) would overflow int64 or
// exceed the 1e15 magnitude budget enforced when clearing coefficient
// denominators. Zero operands are always safe, since their LCM is 0.
func willLCMOverflow(a, b int64) bool {
	if a == 0 || b == 0 {
		return false
	}

	if a > 0 && b > 0 && a > math.MaxInt64/b {
		return true
	}
	if a < 0 && b < 0 && a < math.MaxInt64/b {
		return true
	}
	if (a > 0 && b < 0 && -b > math.MaxInt64/a) || (a < 0 && b > 0 && -a > math.MaxInt64/b) {
		return true
	}

	gcd := gcdInt64(a, b)
	if gcd == 0 {
		return false
	}

	absA := a
	if a < 0 {
		absA = -a
	}
	absB := b
	if b < 0 {
		absB = -b
	}

	if absA/gcd > 1e15/absB {
		return true
	}

	return false
}

// lcmInt64 returns the absolute least common multiple of a and b, or 0 if
// either operand is 0. The caller is responsible for overflow checks; use
// [FindLCMSliceInt64] when working with slices.
func lcmInt64(a, b int64) int64 {
	if a == 0 || b == 0 {
		return 0
	}

	gcd := gcdInt64(a, b)
	result := (a / gcd) * b
	if result < 0 {
		result = -result
	}
	return result
}

// gcdInt64 returns the greatest common divisor of the absolute values of a
// and b, computed with Euclid's algorithm. gcdInt64(0, 0) is 0.
func gcdInt64(a, b int64) int64 {
	if a < 0 {
		a = -a
	}
	if b < 0 {
		b = -b
	}

	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// FindGCDSliceInt64 returns the greatest common divisor of the absolute
// values of nums; the GCD of an empty slice is 1, and 0 is returned only if
// every element is 0. The scan stops as soon as the running GCD reaches 1,
// since it cannot decrease any further.
func FindGCDSliceInt64(nums []int64) int64 {
	if len(nums) == 0 {
		return 1
	}
	result := nums[0]
	if result < 0 {
		result = -result
	}
	for i := 1; i < len(nums); i++ {
		result = gcdInt64(result, nums[i])
		if result == 1 {
			break
		}
	}
	return result
}

// SymmetricDifference returns the elements that occur in only one of the two
// slices: first the elements of slice1 that are missing from slice2, then the
// elements of slice2 that are missing from slice1, each in its original
// order. It detects elements that appear on only one side of a reaction,
// i.e. the set(reactants) ^ set(products). The result is nil when both slices
// contain the same set of elements.
func SymmetricDifference(slice1, slice2 []string) []string {
	set1 := make(map[string]bool)
	set2 := make(map[string]bool)
	for _, v := range slice1 {
		set1[v] = true
	}
	for _, v := range slice2 {
		set2[v] = true
	}

	var result []string
	for _, v := range slice1 {
		if !set2[v] {
			result = append(result, v)
		}
	}
	for _, v := range slice2 {
		if !set1[v] {
			result = append(result, v)
		}
	}

	return result
}

// ReplaceNthOccurrence returns a copy of s in which the n-th (1-based)
// occurrence of old is replaced by new. It returns s unchanged when n <= 0,
// when old is empty or when s contains fewer than n occurrences. It is used
// to insert the "=" between the last reactant and the first product of a
// generated final reaction string.
func ReplaceNthOccurrence(s, old, new string, n int) string {
	if n <= 0 || old == "" {
		return s
	}
	start := 0
	for i := 1; i <= n; i++ {
		index := strings.Index(s[start:], old)
		if index == -1 {
			return s
		}

		if i == n {
			actualIndex := start + index
			return s[:actualIndex] + new + s[actualIndex+len(old):]
		}
		start = start + index + len(old)
	}

	return s
}
