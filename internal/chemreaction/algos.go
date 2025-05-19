package chemreaction

import (
	"fmt"
	"math"
	"slices"

	"gonum.org/v1/gonum/mat"
)

type BalancingAlgos struct {
	ReactionMatrix *mat.Dense
	SeparatorPos   int
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

	var svd mat.SVD
	ok := svd.Factorize(matrix, mat.SVDFull)
	if !ok {
		return nil, fmt.Errorf("SVD factorization failed")
	}

	const tol = 1e-100
	singularValues := make([]float64, min(rows, cols))
	svd.Values(singularValues)

	rank := 0
	for _, val := range singularValues {
		if val > tol {
			rank++
		}
	}

	inverse := mat.NewDense(cols, rows, nil)
	b_ := mat.NewDense(rows, rows, nil)
	for i := range rows {
		b_.Set(i, i, 1.0)
	}

	svd.SolveTo(inverse, b_, rank)

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

func (ba *BalancingAlgos) PPInvAlgorithm() []float64 {
	// Extract reactant and product matrices based on SeparatorPos
	rows, cols := ba.ReactionMatrix.Dims()

	// Extract reactant matrix (columns 0 to SeparatorPos)
	reactantMatrix := mat.NewDense(rows, ba.SeparatorPos, nil)
	for i := 0; i < rows; i++ {
		for j := 0; j < ba.SeparatorPos; j++ {
			reactantMatrix.Set(i, j, ba.ReactionMatrix.At(i, j))
		}
	}

	// Extract product matrix (columns SeparatorPos to end)
	productCols := cols - ba.SeparatorPos
	productMatrix := mat.NewDense(rows, productCols, nil)
	for i := 0; i < rows; i++ {
		for j := 0; j < productCols; j++ {
			productMatrix.Set(i, j, ba.ReactionMatrix.At(i, j+ba.SeparatorPos))
		}
	}

	reactantRows, reactantCols := reactantMatrix.Dims()

	// 1. Calculate Moore-Penrose pseudoinverse of reactant matrix
	mpInverse, err := computePseudoinverse(reactantMatrix)
	if err != nil {
		fmt.Println("Error computing pseudoinverse of reactant matrix:", err)
		return []float64{}
	}

	// 2. Create identity matrix of size reactantRows x reactantRows
	identity := mat.NewDense(reactantRows, reactantRows, nil)
	for i := 0; i < reactantRows; i++ {
		identity.Set(i, i, 1)
	}

	// Calculate (I - A * A^-)
	var reactantMpProduct, tempIdentity mat.Dense
	reactantMpProduct.Mul(reactantMatrix, mpInverse)
	tempIdentity.Sub(identity, &reactantMpProduct)

	// Calculate G = (I - A * A^-) * B
	var gMatrix mat.Dense
	gMatrix.Mul(&tempIdentity, productMatrix)

	// 3. Calculate pseudoinverse of G matrix
	gPinv, err := computePseudoinverse(&gMatrix)
	if err != nil {
		fmt.Println("Error computing pseudoinverse of G matrix:", err)
		return []float64{}
	}

	// Calculate G^- * G
	var gPinvG mat.Dense
	gPinvG.Mul(gPinv, &gMatrix)

	// Get dimensions of G^- * G for creating proper identity matrix
	gPinvGRows, gPinvGCols := gPinvG.Dims()

	// Create identity matrix of same size as G^- * G
	identityGSize := mat.NewDense(gPinvGRows, gPinvGCols, nil)
	for i := 0; i < gPinvGRows; i++ {
		identityGSize.Set(i, i, 1)
	}

	// Calculate (I - G^- * G)
	var yMultiply mat.Dense
	yMultiply.Sub(identityGSize, &gPinvG)

	// Create vector of ones with appropriate size
	_, yMultiplyCols := yMultiply.Dims()
	ones := mat.NewVecDense(yMultiplyCols, nil)
	for i := 0; i < yMultiplyCols; i++ {
		ones.SetVec(i, 1)
	}

	// Calculate y = (I - G^- * G) * ones
	yVector := mat.NewVecDense(yMultiplyCols, nil)
	yVector.MulVec(&yMultiply, ones)

	// 4. Calculate x part

	// Calculate A^- * B
	var mpProduct mat.Dense
	mpProduct.Mul(mpInverse, productMatrix)

	// Calculate A^- * B * y
	var tmpVec mat.VecDense
	tmpVec.MulVec(&mpProduct, yVector)

	// Calculate A^- * A
	var mpA mat.Dense
	mpA.Mul(mpInverse, reactantMatrix)

	// Create identity matrix of size (reactantCols x reactantCols)
	identityA := mat.NewDense(reactantCols, reactantCols, nil)
	for i := 0; i < reactantCols; i++ {
		identityA.Set(i, i, 1)
	}

	// Calculate (I - A^- * A)
	var iMinusMPA mat.Dense
	iMinusMPA.Sub(identityA, &mpA)

	// Create vector of ones with size reactantCols
	vOnes := mat.NewVecDense(reactantCols, nil)
	for i := 0; i < reactantCols; i++ {
		vOnes.SetVec(i, 1)
	}

	// Calculate (I - A^- * A) * v
	var tmpVec2 mat.VecDense
	tmpVec2.MulVec(&iMinusMPA, vOnes)

	// Calculate x = A^-By + (I-A^-A)v
	xVector := mat.NewVecDense(reactantCols, nil)
	xVector.AddVec(&tmpVec, &tmpVec2)

	// Combine x_vector and y_vector into coefficients
	coefs := make([]float64, reactantCols+yMultiplyCols)

	// Copy values from xVector
	for i := 0; i < reactantCols; i++ {
		coefs[i] = xVector.AtVec(i)
	}

	// Copy values from yVector
	for i := 0; i < yMultiplyCols; i++ {
		coefs[reactantCols+i] = yVector.AtVec(i)
	}

	return coefs
}

// computePseudoinverse calculates the Moore-Penrose pseudoinverse of a matrix
// using SVD and SolveTo method
func computePseudoinverse(matrix *mat.Dense) (*mat.Dense, error) {
	rows, cols := matrix.Dims()

	// Create a new SVD factorizer and decompose the matrix
	var svd mat.SVD
	ok := svd.Factorize(matrix, mat.SVDFull)
	if !ok {
		return nil, fmt.Errorf("SVD factorization failed")
	}

	// Get the singular values
	singularValues := make([]float64, min(rows, cols))
	svd.Values(singularValues)

	// Determine numerical rank using tolerance
	const tol = 1e-10 // Tolerance for singular values
	rank := 0
	for _, val := range singularValues {
		if val > tol {
			rank++
		}
	}

	// Create identity matrix of size rows x rows
	b := mat.NewDense(rows, rows, nil)
	for i := 0; i < rows; i++ {
		b.Set(i, i, 1.0)
	}

	// Create the pseudoinverse matrix
	inverse := mat.NewDense(cols, rows, nil)

	// Use SolveTo to compute the pseudoinverse
	svd.SolveTo(inverse, b, rank)

	return inverse, nil
}

// min returns the minimum of two integers
func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func matrixRank(m *mat.Dense) (int, error) {
	rows, cols := m.Dims()
	var svd mat.SVD
	ok := svd.Factorize(m, mat.SVDThin)
	if !ok {
		return -1, fmt.Errorf("SVD factorization failed")
	}
	singularValues := make([]float64, min(rows, cols))
	svd.Values(singularValues)
	rank := 0
	tol := 1e-100
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
