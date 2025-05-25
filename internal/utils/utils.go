package utils

import (
	"math"
	"math/big"
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

type Rational struct {
	Num, Den *big.Int
}

func NewRational(f float64) *Rational {
	if math.IsInf(f, 0) || math.IsNaN(f) {
		return &Rational{big.NewInt(0), big.NewInt(1)}
	}

	rat := big.NewRat(0, 1)
	rat.SetFloat64(f)

	return &Rational{
		Num: new(big.Int).Set(rat.Num()),
		Den: new(big.Int).Set(rat.Denom()),
	}
}

func NewRationalWithLimit(f float64, maxDenominator int64) *Rational {
	if math.IsInf(f, 0) || math.IsNaN(f) {
		return &Rational{big.NewInt(0), big.NewInt(1)}
	}

	if maxDenominator <= 0 {
		maxDenominator = 1000000
	}

	return limitDenominator(f, maxDenominator)
}

func limitDenominator(x float64, maxDen int64) *Rational {
	if x == 0 {
		return &Rational{big.NewInt(0), big.NewInt(1)}
	}

	sign := int64(1)
	if x < 0 {
		sign = -1
		x = -x
	}

	if x == math.Floor(x) {
		return &Rational{big.NewInt(sign * int64(x)), big.NewInt(1)}
	}

	var p0, q0, p1, q1 int64 = 0, 1, 1, 0

	n := x
	for q1 <= maxDen {
		a := int64(math.Floor(n))

		p0, q0, p1, q1 = p1, q1, p1*a+p0, q1*a+q0

		if q1 > maxDen {
			break
		}

		if math.Abs(float64(p1)/float64(q1)-x) < 1e-15 {
			break
		}

		if n == float64(a) {
			break
		}
		n = 1.0 / (n - float64(a))

		if math.IsInf(n, 0) || math.IsNaN(n) {
			break
		}
	}

	if q1 > maxDen {
		p1, q1 = p0, q0
	}

	if q1 == 0 {
		q1 = 1
	}

	return &Rational{big.NewInt(sign * p1), big.NewInt(q1)}
}

func (r *Rational) Simplify() {
	gcd := new(big.Int).GCD(nil, nil, r.Num, r.Den)
	r.Num.Div(r.Num, gcd)
	r.Den.Div(r.Den, gcd)

	if r.Den.Sign() < 0 {
		r.Num.Neg(r.Num)
		r.Den.Neg(r.Den)
	}
}

func gcdBig(a, b *big.Int) *big.Int {
	return new(big.Int).GCD(nil, nil, a, b)
}

func lcmBig(a, b *big.Int) *big.Int {
	if a.Sign() == 0 || b.Sign() == 0 {
		return big.NewInt(0)
	}

	gcd := gcdBig(a, b)
	result := new(big.Int).Mul(a, b)
	result.Div(result, gcd)
	return new(big.Int).Abs(result)
}

func FindLCMSlice(nums []*big.Int) *big.Int {
	if len(nums) == 0 {
		return big.NewInt(1)
	}

	result := new(big.Int).Set(nums[0])
	for i := 1; i < len(nums); i++ {
		result = lcmBig(result, nums[i])
	}
	return result
}

func FindGCDSlice(nums []*big.Int) *big.Int {
	if len(nums) == 0 {
		return big.NewInt(1)
	}

	result := new(big.Int).Set(nums[0])
	for i := 1; i < len(nums); i++ {
		result = gcdBig(result, nums[i])
		if result.Cmp(big.NewInt(1)) == 0 {
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
