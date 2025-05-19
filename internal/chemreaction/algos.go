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

func (ba *BalancingAlgos) PPInvAlgorithm() ([]float64, error) {
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

	// 1. Moore-Penrose pseudoinverse of reactant matrix
	var mpInverse mat.Dense

	// Using SVD method to compute pseudoinverse
	var svd mat.SVD
	ok := svd.Factorize(reactantMatrix, mat.SVDThin)
	if !ok {
		// Handle failed factorization
		return []float64{}, nil
	}

	// Retrieve the components of SVD
	var u, vt mat.Dense
	svd.UTo(&u)
	svd.VTo(&vt)

	// Get singular values
	s := svd.Values(nil)

	// Create diagonal matrix with reciprocals of singular values
	sInv := mat.NewDense(len(s), len(s), nil)
	for i := 0; i < len(s); i++ {
		if s[i] > 1e-15 { // Threshold for numerical stability
			sInv.Set(i, i, 1/s[i])
		} else {
			sInv.Set(i, i, 0)
		}
	}

	// Compute pseudoinverse A^- = V * S^-1 * U^T
	var uT, temp mat.Dense
	uT.CloneFrom(u.T())
	temp.Mul(sInv, &uT)
	mpInverse.Mul(vt.T(), &temp)

	// 2. Create identity matrix of size reactantRows x reactantRows
	identity := mat.NewDense(reactantRows, reactantRows, nil)
	for i := 0; i < reactantRows; i++ {
		identity.Set(i, i, 1)
	}

	// Calculate (I - A * A^-)
	var reactantMpProduct, tempIdentity mat.Dense
	// Note: A * A^- will be of shape (reactantRows x reactantRows)
	reactantMpProduct.Mul(reactantMatrix, &mpInverse)
	tempIdentity.Sub(identity, &reactantMpProduct)

	// Calculate G = (I - A * A^-) * B
	var gMatrix mat.Dense
	gMatrix.Mul(&tempIdentity, productMatrix)

	// 3. Calculate pseudoinverse of G
	var gPinv mat.Dense

	// Use SVD again for G's pseudoinverse
	ok = svd.Factorize(&gMatrix, mat.SVDThin)
	if !ok {
		// Handle failed factorization
		return []float64{}, nil
	}

	// Retrieve SVD components for G
	svd.UTo(&u)
	svd.VTo(&vt)
	s = svd.Values(nil)

	// Create diagonal matrix with reciprocals of G's singular values
	gSInv := mat.NewDense(len(s), len(s), nil)
	for i := 0; i < len(s); i++ {
		if s[i] > 1e-15 { // Threshold for numerical stability
			gSInv.Set(i, i, 1/s[i])
		} else {
			gSInv.Set(i, i, 0)
		}
	}

	// Compute G's pseudoinverse
	uT.CloneFrom(u.T())
	temp.Mul(gSInv, &uT)
	gPinv.Mul(vt.T(), &temp)

	// Get dimensions of G matrix and G pseudoinverse for debugging
	gRows, _ := gMatrix.Dims()
	_, gPinvCols := gPinv.Dims()

	// Debugging: Check dimensions to ensure proper multiplication
	if gPinvCols != gRows {
		// Dimensions don't match, adjust the matrix shape or use a different approach
		return []float64{}, nil
	}

	// Calculate G^- * G
	var gPinvG mat.Dense
	gPinvG.Mul(&gPinv, &gMatrix)

	// Get dimensions of G^- * G
	gPinvGRows, gPinvGCols := gPinvG.Dims()

	// Create identity matrix of same size as G^- * G
	identityGSize := mat.NewDense(gPinvGRows, gPinvGCols, nil)
	for i := 0; i < gPinvGRows; i++ {
		identityGSize.Set(i, i, 1)
	}

	// Calculate (I - G^- * G)
	var yMultiply mat.Dense
	yMultiply.Sub(identityGSize, &gPinvG)

	// Create vector of ones with appropriate size based on yMultiply columns
	_, yMultiplyCols := yMultiply.Dims()
	ones := mat.NewVecDense(yMultiplyCols, nil)
	for i := 0; i < yMultiplyCols; i++ {
		ones.SetVec(i, 1)
	}

	// Calculate y = (I - G^- * G) * ones
	yVector := mat.NewVecDense(yMultiplyCols, nil)
	yVector.MulVec(&yMultiply, ones)

	// 4. Calculate x part

	// Calculate A^- * B (will be reactantCols x productCols)
	var mpProduct mat.Dense
	mpProduct.Mul(&mpInverse, productMatrix)

	// Calculate A^- * B * y (product matrix times y vector)
	var tmpVec mat.VecDense
	tmpVec.MulVec(&mpProduct, yVector)

	// Calculate A^- * A
	var mpA mat.Dense
	mpA.Mul(&mpInverse, reactantMatrix)

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

	return coefs, nil
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
