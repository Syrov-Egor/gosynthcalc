package chemreaction

import (
	"fmt"
	"math"
	"runtime"
	"slices"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

type BalancingAlgos struct {
	ReactionMatrix *mat.Dense
	SeparatorPos   int
	ReactantMatrix *mat.Dense
	ProductMatrix  *mat.Dense
}

func NewBalancingAlgos(reactionMatrix *mat.Dense, separatorPos int) *BalancingAlgos {
	rows, cols := reactionMatrix.Dims()
	reactantMatrix := mat.NewDense(rows, separatorPos, nil)
	for i := range rows {
		for j := range separatorPos {
			reactantMatrix.Set(i, j, reactionMatrix.At(i, j))
		}
	}

	productCols := cols - separatorPos
	productMatrix := mat.NewDense(rows, productCols, nil)
	for i := range rows {
		for j := range productCols {
			productMatrix.Set(i, j, reactionMatrix.At(i, j+separatorPos))
		}
	}
	return &BalancingAlgos{
		ReactionMatrix: reactionMatrix,
		SeparatorPos:   separatorPos,
		ReactantMatrix: reactantMatrix,
		ProductMatrix:  productMatrix,
	}

}

func (b *BalancingAlgos) InvAlgorithm() ([]float64, error) {
	rows, cols := b.ReactionMatrix.Dims()
	reactionMatrix := mat.DenseCopyOf(b.ReactionMatrix)
	var zerosAdded int

	if rows > cols {
		zerosAdded = rows - cols
		newMatrix := mat.NewDense(rows, rows, nil)
		newMatrix.Slice(0, rows, 0, cols).(*mat.Dense).Copy(reactionMatrix)
		reactionMatrix = newMatrix
		cols = rows
	}

	if rows == cols {
		var svd mat.SVD
		ok := svd.Factorize(reactionMatrix, mat.SVDFull)
		if !ok {
			return nil, fmt.Errorf("SVD factorization failed")
		}
		var v mat.Dense
		svd.VTo(&v)
		reactionMatrix = mat.DenseCopyOf(v.T())
	}

	rank, err := matrixRank(reactionMatrix)
	if err != nil {
		return nil, err
	}
	nullity := cols - rank

	augument := mat.NewDense(nullity, cols, nil)
	for i := range nullity {
		augument.Set(i, cols-i-1, 1.0)
	}

	augumentedMatrix := mat.NewDense(rows+nullity, cols, nil)
	augumentedMatrix.Slice(0, rows, 0, cols).(*mat.Dense).Copy(reactionMatrix)
	augumentedMatrix.Slice(rows, rows+nullity, 0, cols).(*mat.Dense).Copy(augument)

	nonZeroRows := findNonZeroRows(augumentedMatrix)
	if len(nonZeroRows) < rows+nullity {
		cleanMatrix := mat.NewDense(len(nonZeroRows), cols, nil)
		for i, rowIdx := range nonZeroRows {
			row := mat.Row(nil, rowIdx, augumentedMatrix)
			cleanMatrix.SetRow(i, row)
		}
		augumentedMatrix = cleanMatrix
	}

	var inversedMatrix mat.Dense

	err = inversedMatrix.Inverse(augumentedMatrix)
	if err != nil {
		return nil, fmt.Errorf("%s", "Matrix inversion failed: "+err.Error())
	}

	r, c := inversedMatrix.Dims()
	vector := make([]float64, r)
	for i := range r {
		vector[i] = inversedMatrix.At(i, c-zerosAdded-1)
	}

	var nonZeroVector []float64
	for _, val := range vector {
		absVal := math.Abs(val)
		if absVal > 1e-10 {
			nonZeroVector = append(nonZeroVector, absVal)
		}
	}

	minVal := slices.Min(nonZeroVector)
	coefs := make([]float64, len(nonZeroVector))
	for i, val := range nonZeroVector {
		coefs[i] = val / minVal
	}

	return coefs, nil
}

func (b *BalancingAlgos) GPInvAlgorithm() ([]float64, error) {
	rows, cols := b.ReactionMatrix.Dims()

	matrix := mat.NewDense(rows, cols, nil)
	matrix.Copy(b.ReactionMatrix)

	for i := range rows {
		for j := b.SeparatorPos; j < cols; j++ {
			matrix.Set(i, j, -matrix.At(i, j))
		}
	}

	inverse, err := computePseudoinverse(matrix)
	if err != nil {
		return nil, err
	}

	identityMatrix := mat.NewDense(cols, cols, nil)
	for i := range cols {
		identityMatrix.Set(i, i, 1.0)
	}

	a := mat.NewVecDense(cols, nil)
	for i := range cols {
		a.SetVec(i, 1.0)
	}

	temp := mat.NewDense(cols, cols, nil)
	temp.Product(inverse, matrix)
	temp.Scale(-1, temp)
	temp.Add(identityMatrix, temp)

	coefs := mat.NewVecDense(cols, nil)
	coefs.MulVec(temp, a)

	result := make([]float64, cols)
	for i := range cols {
		result[i] = coefs.AtVec(i)
	}

	return result, nil
}

func (b *BalancingAlgos) PPInvAlgorithm() ([]float64, error) {
	reactantRows, reactantCols := b.ReactantMatrix.Dims()

	mpInverse, err := computePseudoinverse(b.ReactantMatrix)
	if err != nil {
		return nil, fmt.Errorf("Error computing pseudoinverse of reactant matrix: %s", err)
	}

	identity := mat.NewDense(reactantRows, reactantRows, nil)
	for i := range reactantRows {
		identity.Set(i, i, 1)
	}

	var reactantMpProduct, tempIdentity mat.Dense
	reactantMpProduct.Mul(b.ReactantMatrix, mpInverse)
	tempIdentity.Sub(identity, &reactantMpProduct)

	var gMatrix mat.Dense
	gMatrix.Mul(&tempIdentity, b.ProductMatrix)

	gPinv, err := computePseudoinverse(&gMatrix)
	if err != nil {
		return nil, fmt.Errorf("Error computing pseudoinverse of reactant matrix: %s", err)
	}

	var gPinvG mat.Dense
	gPinvG.Mul(gPinv, &gMatrix)

	gPinvGRows, gPinvGCols := gPinvG.Dims()

	identityGSize := mat.NewDense(gPinvGRows, gPinvGCols, nil)
	for i := range gPinvGRows {
		identityGSize.Set(i, i, 1)
	}

	var yMultiply mat.Dense
	yMultiply.Sub(identityGSize, &gPinvG)

	_, yMultiplyCols := yMultiply.Dims()
	ones := mat.NewVecDense(yMultiplyCols, nil)
	for i := range yMultiplyCols {
		ones.SetVec(i, 1)
	}

	yVector := mat.NewVecDense(yMultiplyCols, nil)
	yVector.MulVec(&yMultiply, ones)

	var mpProduct mat.Dense
	mpProduct.Mul(mpInverse, b.ProductMatrix)

	var tmpVec mat.VecDense
	tmpVec.MulVec(&mpProduct, yVector)

	var mpA mat.Dense
	mpA.Mul(mpInverse, b.ReactantMatrix)

	identityA := mat.NewDense(reactantCols, reactantCols, nil)
	for i := range reactantCols {
		identityA.Set(i, i, 1)
	}

	var iMinusMPA mat.Dense
	iMinusMPA.Sub(identityA, &mpA)

	vOnes := mat.NewVecDense(reactantCols, nil)
	for i := range reactantCols {
		vOnes.SetVec(i, 1)
	}

	var tmpVec2 mat.VecDense
	tmpVec2.MulVec(&iMinusMPA, vOnes)

	xVector := mat.NewVecDense(reactantCols, nil)
	xVector.AddVec(&tmpVec, &tmpVec2)

	coefs := make([]float64, reactantCols+yMultiplyCols)

	for i := range reactantCols {
		coefs[i] = xVector.AtVec(i)
	}

	for i := range yMultiplyCols {
		coefs[reactantCols+i] = yVector.AtVec(i)
	}

	return coefs, nil
}

func (b *BalancingAlgos) Combinatorial(maxCoef uint) []int {
	iMaxCoef := int(maxCoef)
	_, cols := b.ReactionMatrix.Dims()
	gen := NewMultiCombinationGenerator(iMaxCoef, cols)
	arrs := gen.Generate(runtime.GOMAXPROCS(0))

	for arr := range arrs {
		reactantCoefs := arr[:b.SeparatorPos]
		fReactCoefs := make([]float64, len(reactantCoefs))
		for i, val := range reactantCoefs {
			fReactCoefs[i] = float64(val)
		}
		reacSum := mulAndSumMat(b.ReactantMatrix, fReactCoefs)

		productCoefs := arr[b.SeparatorPos:]
		fProdCoefs := make([]float64, len(productCoefs))
		for i, val := range productCoefs {
			fProdCoefs[i] = float64(val)
		}
		prodSum := mulAndSumMat(b.ProductMatrix, fProdCoefs)

		if floats.EqualApprox(reacSum, prodSum, 1e-10) {
			ret := append(reactantCoefs, productCoefs...)
			return ret
		}
	}

	return nil
}

func mulAndSumMat(matrix *mat.Dense, vector []float64) []float64 {
	rows, cols := matrix.Dims()
	res := make([]float64, rows)
	for row := range rows {
		var sum float64
		for el := range cols {
			sum += matrix.At(row, el) * vector[el]
		}
		res[row] = sum
	}
	return res
}

func computePseudoinverse(matrix *mat.Dense) (*mat.Dense, error) {
	rows, cols := matrix.Dims()

	rank, err := matrixRank(matrix)
	if err != nil {
		return nil, err
	}

	b := mat.NewDense(rows, rows, nil)
	for i := range rows {
		b.Set(i, i, 1.0)
	}

	inverse := mat.NewDense(cols, rows, nil)
	var svd mat.SVD

	svd.SolveTo(inverse, b, rank)

	return inverse, nil
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func matrixRank(m *mat.Dense) (int, error) {
	rows, cols := m.Dims()
	var svd mat.SVD
	ok := svd.Factorize(m, mat.SVDFull)
	if !ok {
		return -1, fmt.Errorf("SVD factorization failed")
	}
	singularValues := make([]float64, min(rows, cols))
	svd.Values(singularValues)
	rank := 0
	tol := 1e-12
	for _, val := range singularValues {
		if val > tol {
			rank++
		}
	}

	return rank, nil
}

func findNonZeroRows(m *mat.Dense) []int {
	rows, cols := m.Dims()
	nonZeroRows := []int{}

	for i := range rows {
		isZeroRow := true
		for j := range cols {
			if math.Abs(m.At(i, j)) > 1e-10 {
				isZeroRow = false
				break
			}
		}
		if !isZeroRow {
			nonZeroRows = append(nonZeroRows, i)
		}
	}

	return nonZeroRows
}
