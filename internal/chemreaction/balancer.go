package chemreaction

import (
	"fmt"
	"math/big"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

type Balancer struct {
	reactionMatrix *mat.Dense
	separatorPos   int
	intify         bool
	tolerance      float64
	precision      uint
	bAlgos         *BalancingAlgos
	maxDenom       int
}

type MethodResult struct {
	Method string
	Result []float64
}

func NewBalancer(matrix *mat.Dense, separatorPos int, intify bool, precision uint, tolerance ...float64) *Balancer {
	var tol float64
	if tolerance == nil {
		tol = 1e-8
	} else {
		tol = tolerance[0]
	}

	bAlgos := NewBalancingAlgos(matrix, separatorPos, tol)

	return &Balancer{
		reactionMatrix: matrix,
		separatorPos:   separatorPos,
		intify:         intify,
		tolerance:      tol,
		precision:      precision,
		bAlgos:         bAlgos,
		maxDenom:       1_000_000,
	}
}

func (b *Balancer) reduceCoefficients(coeffs []float64) ([]int64, error) {
	return reduceCoefficientsWithDenomLimit(coeffs, int64(b.maxDenom))
}

func reduceCoefficientsWithDenomLimit(coeffs []float64, maxDenominator int64) ([]int64, error) {
	if len(coeffs) == 0 {
		return []int64{}, nil
	}

	rationals := make([]*utils.Rational, len(coeffs))
	for i, coeff := range coeffs {
		rationals[i] = utils.NewRationalWithLimit(coeff, maxDenominator)
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
			return nil, fmt.Errorf("result too large for int64: coefficient %d", i)
		}

		result[i] = reduced.Int64()
	}

	return result, nil
}

func (b *Balancer) intifyCoefs(coefs []float64, limit int) []float64 {

	reduced, err := b.reduceCoefficients(coefs)
	if err != nil {
		fmt.Println(err)
		return coefs
	}

	for _, coeff := range reduced {
		if int(coeff) > limit {
			fmt.Println(coeff)
			return coefs
		}
	}

	result := make([]float64, len(reduced))
	for i, r := range reduced {
		result[i] = float64(r)
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

func (b *Balancer) calculateByMethod(method string, maxCoef ...uint) ([]float64, error) {
	var coefficients []float64
	var err error
	errm := fmt.Errorf("Can't balance reaction by %s method", method)

	switch method {
	case "inv":
		coefficients, err = b.bAlgos.InvAlgorithm()
		if err != nil {
			return nil, errm
		}
	case "gpinv":
		coefficients, err = b.bAlgos.GPInvAlgorithm()
		if err != nil {
			return nil, errm
		}
	case "ppinv":
		coefficients, err = b.bAlgos.PPInvAlgorithm()
		if err != nil {
			return nil, errm
		}
	case "comb":
		coefficients = b.bAlgos.Combinatorial(maxCoef[0])
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

func (b *Balancer) Inv() ([]float64, error) {
	res, err := b.calculateByMethod("inv")
	if err != nil {
		return nil, err
	}
	return res, nil
}

func (b *Balancer) GPinv() ([]float64, error) {
	res, err := b.calculateByMethod("gpinv")
	if err != nil {
		return nil, err
	}
	return res, nil
}

func (b *Balancer) PPinv() ([]float64, error) {
	res, err := b.calculateByMethod("ppinv")
	if err != nil {
		return nil, err
	}
	return res, nil
}

func (b *Balancer) Comb(maxCoef uint) ([]float64, error) {
	res, err := b.calculateByMethod("comb", maxCoef)
	if err != nil {
		return nil, err
	}
	return res, nil
}

func (b *Balancer) Auto() (MethodResult, error) {
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
