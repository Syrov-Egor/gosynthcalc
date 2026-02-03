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

type SimpleFraction struct {
	Num, Den int64
}

func NewSimpleFraction(f float64, maxDenominator int64) SimpleFraction {
	if math.IsInf(f, 0) || math.IsNaN(f) || f == 0 {
		return SimpleFraction{0, 1}
	}

	sign := int64(1)
	if f < 0 {
		sign = -1
		f = -f
	}

	if f == math.Floor(f) && f < float64(maxDenominator) {
		return SimpleFraction{sign * int64(f), 1}
	}

	if f >= float64(maxDenominator) {
		intPart := min(int64(f), maxDenominator)
		return SimpleFraction{sign * intPart, 1}
	}

	lowerNum, lowerDen := int64(0), int64(1)
	upperNum, upperDen := int64(1), int64(0)

	for {
		mediantNum := lowerNum + upperNum
		mediantDen := lowerDen + upperDen
		if mediantDen > maxDenominator {
			break
		}

		mediantValue := float64(mediantNum) / float64(mediantDen)

		if mediantValue < f {
			lowerNum, lowerDen = mediantNum, mediantDen
		} else if mediantValue > f {
			upperNum, upperDen = mediantNum, mediantDen
		} else {
			return SimpleFraction{sign * mediantNum, mediantDen}
		}

		if lowerDen > 0 && upperDen > 0 {
			lowerVal := float64(lowerNum) / float64(lowerDen)
			upperVal := float64(upperNum) / float64(upperDen)
			if upperVal-lowerVal < 1e-15 {
				break
			}
		}
	}

	bestNum, bestDen := lowerNum, lowerDen
	bestError := math.Abs(f - float64(lowerNum)/float64(lowerDen))

	if upperDen > 0 {
		upperError := math.Abs(f - float64(upperNum)/float64(upperDen))
		if upperError < bestError {
			bestNum, bestDen = upperNum, upperDen
			bestError = upperError
		}
	}

	return SimpleFraction{sign * bestNum, bestDen}
}

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
