package utils

import (
	"math"
	"strings"
)

func StringCounter(s string) map[string]int {
	counts := make(map[string]int)
	for _, char := range s {
		counts[string(char)]++
	}
	return counts
}

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

func SumFloatS(s []float64) float64 {
	var sum float64 = 0.0
	for _, el := range s {
		sum += el
	}
	return sum
}

func RoundFloat(val float64, precision uint) float64 {
	ratio := math.Pow(10, float64(precision))
	return math.Round(val*ratio) / ratio
}

func RoundFloatS(s []float64, precision uint) []float64 {
	res := make([]float64, len(s))
	for i, val := range s {
		res[i] = RoundFloat(val, precision)
	}
	return res
}

// Simple fraction struct for the intify process
type SimpleFraction struct {
	Num, Den int64
}

// Create fraction from float64 using limit_denominator approach similar to Python
func NewSimpleFraction(f float64, maxDenominator int64) SimpleFraction {
	if math.IsInf(f, 0) || math.IsNaN(f) || f == 0 {
		return SimpleFraction{0, 1}
	}

	// Handle negative numbers
	sign := int64(1)
	if f < 0 {
		sign = -1
		f = -f
	}

	// If it's already an integer
	if f == math.Floor(f) && f < float64(maxDenominator) {
		return SimpleFraction{sign * int64(f), 1}
	}

	// Use Python's Fraction.limit_denominator() algorithm
	// This is simpler than the continued fraction approach
	bestNum, bestDen := int64(0), int64(1)
	bestError := math.Abs(f)

	for den := int64(1); den <= maxDenominator; den++ {
		num := int64(math.Round(f * float64(den)))
		if num == 0 && f != 0 {
			continue
		}

		error := math.Abs(f - float64(num)/float64(den))
		if error < bestError {
			bestError = error
			bestNum = num
			bestDen = den
		}

		// If we found an exact match, stop
		if error < 1e-15 {
			break
		}
	}

	return SimpleFraction{sign * bestNum, bestDen}
}

// Find LCM of a slice of integers
func FindLCMSliceInt64(nums []int64) int64 {
	if len(nums) == 0 {
		return 1
	}

	result := nums[0]
	for i := 1; i < len(nums); i++ {
		// Check for potential overflow before calculation
		if willLCMOverflow(result, nums[i]) {
			return -1
		}
		result = lcmInt64(result, nums[i])
	}
	return result
}

// Helper function to check if LCM calculation will overflow
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

// Find GCD of a slice of integers
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
