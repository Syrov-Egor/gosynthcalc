package chemreaction

import (
	"fmt"
	"math/big"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
	"gonum.org/v1/gonum/mat"
)

type Number interface {
	int | int64 | float64
}

type Balancer struct {
	reactionMatrix *mat.Dense
	separatorPos   int
	intify         bool
	tolerance      float64
	bAlgos         *BalancingAlgos
}

func NewBalancer(matrix *mat.Dense, separatorPos int, intify bool, tolerance ...float64) *Balancer {
	var tol float64
	if tolerance == nil {
		tol = 1e-12
	} else {
		tol = tolerance[0]
	}

	bAlgos := NewBalancingAlgos(matrix, separatorPos, tol)

	return &Balancer{
		reactionMatrix: matrix,
		separatorPos:   separatorPos,
		intify:         intify,
		tolerance:      tol,
		bAlgos:         bAlgos,
	}
}

func reduceCoefs(coefs []float64) ([]int64, error) {
	rationals := make([]*utils.Rational, len(coefs))
	for i, coef := range coefs {
		rationals[i] = utils.NewRational(coef)
		rationals[i].Simplify()
	}
	denominators := make([]*big.Int, len(rationals))
	for i, rat := range rationals {
		denominators[i] = new(big.Int).Set(rat.Den)
	}

	lcm := utils.FindLCMSlice(denominators)

	intVals := make([]*big.Int, len(rationals))
	for i, rat := range rationals {
		quotient := new(big.Int).Div(lcm, rat.Den)
		intVals[i] = new(big.Int).Mul(rat.Num, quotient)
	}

	gcd := utils.FindGCDSlice(intVals)

	result := make([]int64, len(intVals))
	for i, val := range intVals {
		reduced := new(big.Int).Div(val, gcd)

		if !reduced.IsInt64() {
			return nil, fmt.Errorf("Result too large for int64: coefficient %d", i)
		}

		result[i] = reduced.Int64()
	}

	return result, nil
}

func intifyCoefs[Num Number](coefs []Num, limit int) []Num {
	switch any(coefs).(type) {

	case []float64:
		float64Coefs := make([]float64, len(coefs))
		for i, c := range coefs {
			float64Coefs[i] = float64(c)
		}
		reduced, err := reduceCoefs(float64Coefs)
		if err != nil {
			return coefs
		}
		for _, coeff := range reduced {
			if int(coeff) > limit {
				return coefs
			}
		}
		result := make([]Num, len(reduced))
		for i, r := range reduced {
			result[i] = Num(float64(r))
		}
		return result

	case []int, []int64:
		return coefs

	default:
		return coefs
	}
}
