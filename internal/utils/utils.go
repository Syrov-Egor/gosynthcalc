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
