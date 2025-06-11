package chemreaction

import (
	"fmt"
	"math"
	"runtime"
	"slices"
	"sync"
	"sync/atomic"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

type balancingAlgos struct {
	ReactionMatrix *mat.Dense
	SeparatorPos   int
	ReactantMatrix *mat.Dense
	ProductMatrix  *mat.Dense
	ReactantRows   int
	ProductRows    int
	ReactantCols   int
	ProductCols    int
	Tolerance      float64
}

func newBalancingAlgos(reactionMatrix *mat.Dense, separatorPos int, tolerance ...float64) *balancingAlgos {
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

	var tol float64
	if tolerance == nil {
		tol = 1e-12
	} else {
		tol = tolerance[0]
	}

	return &balancingAlgos{
		ReactionMatrix: reactionMatrix,
		SeparatorPos:   separatorPos,
		ReactantMatrix: reactantMatrix,
		ProductMatrix:  productMatrix,
		ReactantRows:   rows,
		ProductRows:    rows,
		ReactantCols:   separatorPos,
		ProductCols:    productCols,
		Tolerance:      tol,
	}
}

func (b *balancingAlgos) invAlgorithm() ([]float64, error) {
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

	rank, _, err := matrixRank(reactionMatrix, b.Tolerance)
	if err != nil {
		return nil, err
	}
	nullity := cols - rank

	var augumentedMatrix *mat.Dense

	if nullity > 0 {
		augument := mat.NewDense(nullity, cols, nil)
		for i := range nullity {
			augument.Set(i, cols-i-1, 1.0)
		}
		augumentedMatrix = mat.NewDense(rows+nullity, cols, nil)
		augumentedMatrix.Slice(0, rows, 0, cols).(*mat.Dense).Copy(reactionMatrix)
		augumentedMatrix.Slice(rows, rows+nullity, 0, cols).(*mat.Dense).Copy(augument)
	} else {
		augumentedMatrix = mat.NewDense(rows, cols, nil)
		augumentedMatrix.Slice(0, rows, 0, cols).(*mat.Dense).Copy(reactionMatrix)
	}

	nonZeroRows := findNonZeroRows(augumentedMatrix, b.Tolerance)
	if len(nonZeroRows) < rows+nullity {
		cleanMatrix := mat.NewDense(len(nonZeroRows), cols, nil)
		for i, rowIdx := range nonZeroRows {
			row := mat.Row(nil, rowIdx, augumentedMatrix)
			cleanMatrix.SetRow(i, row)
		}
		augumentedMatrix = cleanMatrix
	}

	aRows, aCols := augumentedMatrix.Dims()

	if aRows != aCols {
		return nil, fmt.Errorf("Singular matrix")
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
		if absVal > b.Tolerance {
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

func (b *balancingAlgos) gPInvAlgorithm() ([]float64, error) {
	rows, cols := b.ReactionMatrix.Dims()

	matrix := mat.NewDense(rows, cols, nil)
	matrix.Copy(b.ReactionMatrix)

	for i := range rows {
		for j := b.SeparatorPos; j < cols; j++ {
			matrix.Set(i, j, -matrix.At(i, j))
		}
	}

	inverse, err := computePseudoinverse(matrix, b.Tolerance)
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

func (b *balancingAlgos) pPInvAlgorithm() ([]float64, error) {
	reactantRows, reactantCols := b.ReactantMatrix.Dims()

	mpInverse, err := computePseudoinverse(b.ReactantMatrix, b.Tolerance)
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

	gPinv, err := computePseudoinverse(&gMatrix, b.Tolerance)
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

func (b *balancingAlgos) combinatorial(maxCoef uint) []float64 {
	iMaxCoef := int(maxCoef)
	_, cols := b.ReactionMatrix.Dims()
	numWorkers := runtime.GOMAXPROCS(0)
	gen := newMultiCombinationGenerator(iMaxCoef, cols)
	combinations := gen.generate(numWorkers)

	resultChan := make(chan []int, 1)
	var activeWorkers int32
	var wg sync.WaitGroup

	for range numWorkers {
		wg.Add(1)
		go func() {
			defer wg.Done()
			reacSum := make([]float64, b.ReactantRows)
			prodSum := make([]float64, b.ProductRows)

			for arr := range combinations {
				select {
				case result := <-resultChan:
					resultChan <- result
					return
				default:
				}

				atomic.AddInt32(&activeWorkers, 1)

				reactantCoefs := arr[:b.SeparatorPos]
				productCoefs := arr[b.SeparatorPos:]

				mulAndSum(b.ReactantMatrix, reactantCoefs, reacSum, b.ReactantRows, b.ReactantCols)
				mulAndSum(b.ProductMatrix, productCoefs, prodSum, b.ProductRows, b.ProductCols)

				if floats.EqualApprox(reacSum, prodSum, b.Tolerance) {
					solution := arr
					select {
					case resultChan <- solution:
					default:
					}
				}
				atomic.AddInt32(&activeWorkers, -1)
			}
		}()
	}
	go func() {
		wg.Wait()
		close(resultChan)
	}()
	result, ok := <-resultChan
	if !ok {
		return nil
	}
	fmt.Println()
	toFloat := make([]float64, len(result))
	for i, r := range result {
		toFloat[i] = float64(r)
	}
	return toFloat
}

func mulAndSum(matrix *mat.Dense, vector []int, result []float64, rows int, cols int) {
	for row := range rows {
		result[row] = 0
		for col := range cols {
			result[row] += matrix.At(row, col) * float64(vector[col])
		}
	}
}

func computePseudoinverse(matrix *mat.Dense, tol float64) (*mat.Dense, error) {
	rows, cols := matrix.Dims()

	rank, svd, err := matrixRank(matrix, tol)
	if err != nil {
		return nil, err
	}

	b := mat.NewDense(rows, rows, nil)
	for i := range rows {
		b.Set(i, i, 1.0)
	}

	inverse := mat.NewDense(cols, rows, nil)

	svd.SolveTo(inverse, b, rank)

	return inverse, nil
}

func matrixRank(m *mat.Dense, tol float64) (int, mat.SVD, error) {
	rows, cols := m.Dims()
	var svd mat.SVD
	ok := svd.Factorize(m, mat.SVDFull)
	if !ok {
		return -1, svd, fmt.Errorf("SVD factorization failed")
	}
	singularValues := make([]float64, min(rows, cols))
	svd.Values(singularValues)
	rank := 0
	for _, val := range singularValues {
		if val > tol {
			rank++
		}
	}

	return rank, svd, nil
}

func findNonZeroRows(m *mat.Dense, tol float64) []int {
	rows, cols := m.Dims()
	nonZeroRows := []int{}

	for i := range rows {
		isZeroRow := true
		for j := range cols {
			if math.Abs(m.At(i, j)) > tol {
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

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
