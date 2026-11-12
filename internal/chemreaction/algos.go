package chemreaction

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

// balancingAlgos implements the legacy matrix heuristics used to balance a
// chemical reaction. The reaction matrix is split at SeparatorPos into the
// reactant and product halves that every algorithm works with. The exact
// nullity-one solver lives in nullity.go and does not use these heuristics.
//
// Fields:
//
//   - ReactionMatrix is the whole matrix (elements × compounds);
//   - SeparatorPos is the column index of the first product;
//   - ReactantMatrix and ProductMatrix are the left and right halves;
//   - ReactantRows, ProductRows, ReactantCols and ProductCols are their
//     dimensions (both row counts equal the number of elements);
//   - Tolerance is the legacy float tolerance used by the heuristics; rank
//     truncation treats it as relative, and balance checks scale it by the
//     coefficient magnitude. The exact router does not use it.
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

// newBalancingAlgos splits reactionMatrix into its reactant columns
// (0..separatorPos) and product columns (separatorPos..) and returns the
// algorithm collection for it. The tolerance defaults to 1e-12 when omitted;
// [balancer.newBalancer] always passes the balancer's own tolerance.
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

// gPInvAlgorithm computes the coefficients with the matrix general
// pseudoinverse algorithm proposed by
// [Risteski](http://koreascience.or.kr/article/JAKO201314358624990.page).
// It stacks the reactant matrix and the negated product matrix, computes the
// Moore-Penrose pseudoinverse A⁺ and projects the ones vector onto the null
// space:
//
//	x = (I - A⁺ A) 1.
//
// The projection is always a null vector, but it is not guaranteed to be
// strictly positive. For nullity 1 and a positive null direction it is a
// positive multiple of the solution; for nullity 2 or more it can leave the
// positive orthant. A non-positive candidate is reported as
// [ErrMethodNotPositive]; a failed heuristic does not prove infeasibility.
//
// rank is the exact signed-matrix rank from the balancer's exact analysis.
// Using it instead of an absolute SVD tolerance makes the truncation
// deterministic. The returned vector contains raw float values; rounding,
// intification and validation happen in [balancer.calculateLegacyByMethod].
func (b *balancingAlgos) gPInvAlgorithm(rank int) ([]float64, error) {
	rows, cols := b.ReactionMatrix.Dims()

	matrix := mat.NewDense(rows, cols, nil)
	matrix.Copy(b.ReactionMatrix)

	for i := range rows {
		for j := b.SeparatorPos; j < cols; j++ {
			matrix.Set(i, j, -matrix.At(i, j))
		}
	}

	inverse, err := computePseudoinverse(matrix, rank)
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

	if !allStrictlyPositive(result, b.Tolerance) {
		return nil, fmt.Errorf("%w: gpinv projection is not strictly positive", ErrMethodNotPositive)
	}
	return result, nil
}

// pPInvAlgorithm computes the coefficients with the matrix partial
// pseudoinverse algorithm also proposed by
// [Risteski](https://www.koreascience.or.kr/article/JAKO200802727293429.page).
// It uses the general solution
//
//	G = (I - A A⁺) B,   y = (I - G⁺ G) u,   x = A⁺ B y + (I - A⁺ A) v
//
// with the arbitrary choice u = v = 1, and returns [x, y].
//
// The choice u = v = 1 does not guarantee positivity, so the method is a
// heuristic. When G is zero, G⁺ is zero and (I - G⁺G) is the identity, so
// y = u directly; this branch is handled explicitly instead of asking for the
// pseudoinverse of a zero matrix. A non-positive candidate is reported as
// [ErrMethodNotPositive], never as infeasibility.
//
// reactantRank is the exact rank of the reactant matrix. The rank of the
// intermediate G matrix is estimated with the legacy relative-tolerance SVD
// helper because G is already a floating-point intermediate.
func (b *balancingAlgos) pPInvAlgorithm(reactantRank int) ([]float64, error) {
	reactantRows, reactantCols := b.ReactantMatrix.Dims()

	mpInverse, err := computePseudoinverse(b.ReactantMatrix, reactantRank)
	if err != nil {
		return nil, fmt.Errorf("error computing pseudoinverse of reactant matrix: %w", err)
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

	// G = (I - A A⁺)B is mathematically zero when B lies in the column
	// space of A. In floating point it can be tiny non-zero noise, so use a
	// scale-aware zero test before falling back to the relative SVD rank.
	gRank := 0
	if !denseIsZero(&gMatrix, b.Tolerance*max(1.0, denseMaxAbs(b.ProductMatrix))) {
		gRank, _, err = matrixRank(&gMatrix, b.Tolerance)
		if err != nil {
			return nil, fmt.Errorf("error computing rank of G: %w", err)
		}
	}

	var gPinv *mat.Dense
	if gRank == 0 {
		// G = 0  =>  G⁺ = 0  =>  (I - G⁺G) = I  =>  y = u.
		gPinv = mat.NewDense(b.ProductCols, reactantRows, nil)
	} else {
		gPinv, err = computePseudoinverse(&gMatrix, gRank)
		if err != nil {
			return nil, fmt.Errorf("error computing pseudoinverse of G: %w", err)
		}
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

	if !allStrictlyPositive(coefs, b.Tolerance) {
		return nil, fmt.Errorf("%w: ppinv result is not strictly positive", ErrMethodNotPositive)
	}
	return coefs, nil
}

// computePseudoinverse returns the Moore-Penrose pseudoinverse of matrix as
// a cols×rows dense matrix, truncated to the supplied rank. A rank of zero
// returns the zero matrix: this is the mathematically correct pseudoinverse
// and makes zero matrices usable by the legacy heuristics.
func computePseudoinverse(matrix *mat.Dense, rank int) (*mat.Dense, error) {
	if matrix == nil {
		return nil, fmt.Errorf("%w: nil matrix", ErrInvalidMatrix)
	}
	rows, cols := matrix.Dims()
	maxRank := min(rows, cols)
	if rank < 0 || rank > maxRank {
		return nil, fmt.Errorf("rank %d out of range [0,%d]", rank, maxRank)
	}
	if rank == 0 {
		return mat.NewDense(cols, rows, nil), nil
	}

	var svd mat.SVD
	if !svd.Factorize(matrix, mat.SVDFull) {
		return nil, fmt.Errorf("SVD factorization failed")
	}

	b := mat.NewDense(rows, rows, nil)
	for i := range rows {
		b.Set(i, i, 1.0)
	}

	inverse := mat.NewDense(cols, rows, nil)
	svd.SolveTo(inverse, b, rank)
	return inverse, nil
}

// matrixRank returns the rank of m — the number of singular values above a
// relative threshold — together with the SVD factorization, so callers that
// also need the decomposition can reuse it. The threshold is
// tol * max(rows,cols) * σ_max, which is scale-aware unlike the old absolute
// tolerance. It is used only by legacy diagnostics; the exact router never
// calls it.
func matrixRank(m *mat.Dense, tol float64) (int, mat.SVD, error) {
	if m == nil {
		return -1, mat.SVD{}, fmt.Errorf("%w: nil matrix", ErrInvalidMatrix)
	}
	rows, cols := m.Dims()
	var svd mat.SVD
	ok := svd.Factorize(m, mat.SVDFull)
	if !ok {
		return -1, svd, fmt.Errorf("SVD factorization failed")
	}
	singularValues := make([]float64, min(rows, cols))
	svd.Values(singularValues)
	if len(singularValues) == 0 {
		return 0, svd, nil
	}

	if tol < 0 {
		tol = 0
	}
	threshold := tol * float64(max(rows, cols)) * singularValues[0]
	rank := 0
	for _, val := range singularValues {
		if val > 0 && val > threshold {
			rank++
		}
	}
	return rank, svd, nil
}

// denseMaxAbs returns the largest absolute entry of m (0 for an empty
// matrix).
func denseMaxAbs(m *mat.Dense) float64 {
	if m == nil {
		return 0
	}
	rows, cols := m.Dims()
	maxAbs := 0.0
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			if a := math.Abs(m.At(i, j)); a > maxAbs {
				maxAbs = a
			}
		}
	}
	return maxAbs
}

// denseIsZero reports whether every entry of m is at most tol in absolute
// value. NaN and infinities make it return false.
func denseIsZero(m *mat.Dense, tol float64) bool {
	if m == nil {
		return true
	}
	if tol < 0 {
		tol = 0
	}
	rows, cols := m.Dims()
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			v := m.At(i, j)
			if math.IsNaN(v) || math.IsInf(v, 0) || math.Abs(v) > tol {
				return false
			}
		}
	}
	return true
}

// allStrictlyPositive reports whether every entry is positive relative to a
// scaled tolerance. NaN and infinities are rejected.
func allStrictlyPositive(v []float64, tol float64) bool {
	if len(v) == 0 {
		return false
	}
	scale := 1.0
	for _, x := range v {
		if math.IsNaN(x) || math.IsInf(x, 0) {
			return false
		}
		scale += math.Abs(x)
	}
	if tol < 0 {
		tol = 0
	}
	atol := tol * scale
	for _, x := range v {
		if x <= atol {
			return false
		}
	}
	return true
}
