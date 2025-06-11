package chemreaction

import (
	"fmt"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

type balancer struct {
	reactionMatrix *mat.Dense
	separatorPos   int
	intify         bool
	tolerance      float64
	precision      uint
	bAlgos         *balancingAlgos
	maxDenom       int
}

type MethodResult struct {
	Method string
	Result []float64
}

func newBalancer(matrix *mat.Dense, separatorPos int, intify bool, precision uint, tolerance ...float64) *balancer {
	var tol float64
	if tolerance == nil {
		tol = 1e-8
	} else {
		tol = tolerance[0]
	}

	bAlgos := newBalancingAlgos(matrix, separatorPos, tol)

	return &balancer{
		reactionMatrix: matrix,
		separatorPos:   separatorPos,
		intify:         intify,
		tolerance:      tol,
		precision:      precision,
		bAlgos:         bAlgos,
		maxDenom:       1_000_000,
	}
}

func (b *balancer) intifyCoefs(coefs []float64, limit int) []float64 {
	initialCoefficients := make([]float64, len(coefs))
	copy(initialCoefficients, coefs)

	fractions := make([]utils.SimpleFraction, len(coefs))
	denominators := make([]int64, len(coefs))

	for i, coef := range coefs {
		fractions[i] = utils.NewSimpleFraction(coef, int64(b.maxDenom))
		denominators[i] = fractions[i].Den
	}

	lcm := utils.FindLCMSliceInt64(denominators)
	if lcm < 0 || lcm > 1e15 {
		return initialCoefficients
	}

	vals := make([]int64, len(fractions))
	for i, frac := range fractions {
		if frac.Den == 0 {
			return initialCoefficients
		}
		vals[i] = frac.Num * (lcm / frac.Den)

		if vals[i] < 0 && frac.Num > 0 {
			return initialCoefficients
		}
	}

	gcd := utils.FindGCDSliceInt64(vals)
	if gcd == 0 {
		return initialCoefficients
	}

	coefficients := make([]int64, len(vals))
	for i, val := range vals {
		coefficients[i] = val / gcd
	}

	for _, coeff := range coefficients {
		if coeff < 0 {
			coeff = -coeff
		}
		if int(coeff) > limit {
			return initialCoefficients
		}
	}

	result := make([]float64, len(coefficients))
	for i, coeff := range coefficients {
		result[i] = float64(coeff)
	}

	return result
}

func isReactionBalanced(reactantMatrix *mat.Dense, productMatrix *mat.Dense, coefs []float64, atol float64) bool {
	reactantRows, reactantCols := reactantMatrix.Dims()
	productRows, productCols := productMatrix.Dims()
	separatorPos := reactantCols
	reactantCoefs := coefs[:separatorPos]
	productCoefs := coefs[separatorPos:]
	reacSum := make([]float64, reactantRows)
	prodSum := make([]float64, productRows)
	mulAndSumFl(reactantMatrix, reactantCoefs, reacSum, reactantRows, reactantCols)
	mulAndSumFl(productMatrix, productCoefs, prodSum, productRows, productCols)

	if floats.EqualApprox(reacSum, prodSum, atol) {
		return true
	}
	return false
}

func (b *balancer) calculateByMethod(method string, maxCoef ...uint) ([]float64, error) {
	var coefficients []float64
	var err error
	errm := fmt.Errorf("Can't balance reaction by %s method", method)

	switch method {
	case "inv":
		coefficients, err = b.bAlgos.invAlgorithm()
		if err != nil {
			return nil, errm
		}
	case "gpinv":
		coefficients, err = b.bAlgos.gPInvAlgorithm()
		if err != nil {
			return nil, errm
		}
	case "ppinv":
		coefficients, err = b.bAlgos.pPInvAlgorithm()
		if err != nil {
			return nil, errm
		}
	case "comb":
		coefficients = b.bAlgos.combinatorial(maxCoef[0])
		if coefficients == nil {
			return nil, errm
		}
	default:
		return nil, fmt.Errorf("No method %s", method)
	}
	coefficients = utils.RoundFloatS(coefficients, b.precision+2)
	_, matrLength := b.reactionMatrix.Dims()

	if len(coefficients) == matrLength &&
		allPositive(coefficients) &&
		isReactionBalanced(
			b.bAlgos.ReactantMatrix,
			b.bAlgos.ProductMatrix,
			coefficients,
			b.tolerance,
		) {
		if b.intify {
			coefficients = b.intifyCoefs(coefficients, b.maxDenom)
		}
		return coefficients, nil
	}

	return nil, fmt.Errorf("Wrong coefficients")
}

func (b *balancer) Inv() ([]float64, error) {
	res, err := b.calculateByMethod("inv")
	if err != nil {
		return nil, err
	}
	return res, nil
}

func (b *balancer) GPinv() ([]float64, error) {
	res, err := b.calculateByMethod("gpinv")
	if err != nil {
		return nil, err
	}
	return res, nil
}

func (b *balancer) PPinv() ([]float64, error) {
	res, err := b.calculateByMethod("ppinv")
	if err != nil {
		return nil, err
	}
	return res, nil
}

func (b *balancer) Comb(maxCoef uint) ([]float64, error) {
	res, err := b.calculateByMethod("comb", maxCoef)
	if err != nil {
		return nil, err
	}
	return res, nil
}

func (b *balancer) Auto() (MethodResult, error) {
	var coefs []float64
	var err error

	coefs, err = b.Inv()
	if err == nil {
		return MethodResult{Method: "inverse", Result: coefs}, nil
	}
	coefs, err = b.GPinv()
	if err == nil {
		return MethodResult{Method: "general pseudoinverse", Result: coefs}, nil
	}
	coefs, err = b.PPinv()
	if err == nil {
		return MethodResult{Method: "partial pseudoinverse", Result: coefs}, nil
	}

	return MethodResult{Method: "", Result: nil},
		fmt.Errorf("Can't balance this reaction by any method")
}

func mulAndSumFl(matrix *mat.Dense, vector []float64, result []float64, rows int, cols int) {
	for row := range rows {
		result[row] = 0
		for col := range cols {
			result[row] += matrix.At(row, col) * vector[col]
		}
	}
}

func allPositive(coefs []float64) bool {
	for _, coef := range coefs {
		if coef <= 0 {
			return false
		}
	}
	return true
}
