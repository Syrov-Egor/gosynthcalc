package chemreaction

import (
	"context"
	"fmt"
	"math"
	"runtime"
	"slices"
	"sync"
	"sync/atomic"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

// balancingAlgos implements the four algorithms used to balance a chemical
// reaction. The reaction matrix is split at SeparatorPos into the reactant and product
// halves that every algorithm works with.
//
// Fields:
//
//   - ReactionMatrix is the whole matrix (elements × compounds);
//   - SeparatorPos is the column index of the first product;
//   - ReactantMatrix and ProductMatrix are the left and right halves;
//   - ReactantRows, ProductRows, ReactantCols and ProductCols are their
//     dimensions (both row counts equal the number of elements);
//   - Tolerance is the absolute tolerance used for ranks, pseudoinverses and
//     balance comparisons.
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

// invAlgorithm computes the coefficients with the reaction matrix inverse
// algorithm proposed by [Thorne](https://arxiv.org/abs/1110.4321). The
// calculation is based on the nullity, or dimensionality, of the matrix.
//
// The algorithm can be described in steps:
//
//  1. If the number of rows is greater than the number of columns, add zero
//     columns until the matrix becomes square (Note: this is a modification
//     of the original Thorne method described in the article).
//
//  2. If the reaction matrix is square (which means that the number of atoms
//     involved is equal to the number of compounds) then turn the matrix
//     into its row-echelon form by singular value decomposition.
//
//  3. Calculate the nullity of the matrix, which is basically the number of
//     compounds minus the rank of the matrix.
//
//  4. Create a matrix augmented by a nullity number of rows of a flipped
//     identity matrix. If any rows are zeros, replace them with identity
//     matrix rows.
//
//  5. Invert the augmented matrix.
//
//  6. Extract the column of the inverse that corresponds to the last
//     original column (the rightmost one, before the zero columns added in
//     step 1) and transpose it.
//
//  7. Normalize this value with the absolute min value of the vector.
//
//  8. Round up float operations errors (done by
//     [balancer.calculateByMethod] once the algorithm returns).
//
// The absolute values of this vector are the coefficients of the reaction.
//
// Note: while this method works great for reactions with 0 and 1 nullity, it
// generally cannot work with nullities 2 and higher. Thorne claims that for
// higher nullities a nullity number of vectors should be extracted, and each
// of them contains a set of correct coefficients. However, if the number of
// rows in the flipped augmentation identity matrix is 2 or more, one can
// easily see that each vector will contain nullity-1 zeroes, therefore they
// cannot be a correct vector of coefficients.
//
// Go implementation counts the singular values above Tolerance. It returns
// an error when the factorization or the inversion of
// the augmented matrix fails.
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
		return nil, fmt.Errorf("singular matrix")
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

// gPInvAlgorithm computes the coefficients with the matrix general
// pseudoinverse algorithm proposed by
// [Risteski](http://koreascience.or.kr/article/JAKO201314358624990.page).
// There are other articles and methods of chemical equation balancing by
// this author, however, this particular algorithm seems to be the most
// convenient for matrix calculations. The algorithm can be described in
// steps:
//
//  1. Stack reactant matrix and negative product matrix.
//
//  2. Calculate MP pseudoinverse of this matrix.
//
//  3. Calculate coefficients by formula:
//     x = (I – A⁺A)a, where x is the coefficients vector,
//     I - identity matrix, A⁺ - MP inverse, A - matrix,
//     a - arbitrary vector (in this case, vector of ones).
//
// Note: this method is more general than Thorne's method, although it has
// some peculiarities of its own. First of all, the output of this method is
// a float array, so, to generate an int coefs list, it needs to be converted,
// which does not always lead to a good result. Secondly, MP pseudoinverse is
// sensitive to row order in the reaction matrix. The rows should be ordered
// by atoms appearance in the reaction string.
//
// The returned vector still contains the raw (possibly negative) values;
// rounding, intification and validation happen in [balancer.calculateByMethod].
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

// pPInvAlgorithm computes the coefficients with the matrix partial
// pseudoinverse algorithm also proposed by
// [Risteski](https://www.koreascience.or.kr/article/JAKO200802727293429.page).
// The method is founded on the virtue of the solution of a Diophantine
// matrix equation by using the Moore-Penrose pseudoinverse matrix.
//
// The algorithm can be described in steps:
//
//  1. Take the Moore-Penrose pseudoinverse of the reactant matrix.
//
//  2. Create a G matrix in the form of (I-AA^-)B, where I is the identity
//     matrix, A is the reactant matrix, A^- is the MP pseudoinverse of A and
//     B is the product matrix.
//
//  3. Then, the vector y (coefficients of products) is equal to (I-G^-G)u.
//
//  4. Vector x (coefficients of reactants) is equal to A^-By + (I-A^-A)v,
//     where u and v are columns of ones.
//
// The result is [x, y]: the reactant coefficients followed by the product
// coefficients.
//
// Note: while this algorithm and [balancingAlgos.gPInvAlgorithm] are very
// similar, there are some differences in output results. This method exists
// mostly for legacy purposes, like balancing some reactions according to
// Risteski.
func (b *balancingAlgos) pPInvAlgorithm() ([]float64, error) {
	reactantRows, reactantCols := b.ReactantMatrix.Dims()

	mpInverse, err := computePseudoinverse(b.ReactantMatrix, b.Tolerance)
	if err != nil {
		return nil, fmt.Errorf("error computing pseudoinverse of reactant matrix: %s", err)
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
		return nil, fmt.Errorf("error computing pseudoinverse of reactant matrix: %s", err)
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

// combinatorial finds a solution of a Diophantine matrix equation by simply
// enumerating candidate coefficient vectors: one coefficient per compound,
// each between 1 and maxCoef, is generated by
// [multiCombinationGenerator] and tested until the weighted reactant and
// product column sums agree within Tolerance.
//
// The candidates are produced and checked by GOMAXPROCS workers in parallel;
// the first solution found wins, so smaller coefficients are likely to be
// reported first. The search honours ctx: it stops when the context is
// cancelled and returns nil in that case, as well as when no combination
// balances the reaction.
//
// The method is only practical for reactions with a small number of
// compounds, since the search space grows as maxCoef**number_of_compounds.
// Progress is printed to stdout by the generator.
func (b *balancingAlgos) combinatorial(ctx context.Context, maxCoef uint) []float64 {
	iMaxCoef := int(maxCoef)
	_, cols := b.ReactionMatrix.Dims()
	numWorkers := runtime.GOMAXPROCS(0)
	gen := newMultiCombinationGenerator(iMaxCoef, cols)
	combinations := gen.generate(ctx, numWorkers)

	resultChan := make(chan []int, 1)
	var activeWorkers atomic.Int32
	var wg sync.WaitGroup

	for range numWorkers {
		wg.Go(func() {
			reacSum := make([]float64, b.ReactantRows)
			prodSum := make([]float64, b.ProductRows)

			for {
				select {
				case <-ctx.Done(): // Handle cancellation
					return
				case arr, ok := <-combinations:
					if !ok {
						return
					}
					for arr = range combinations {
						select {
						case result := <-resultChan:
							resultChan <- result
							return
						default:
						}

						activeWorkers.Add(1)

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
						activeWorkers.Add(-1)
					}
				}
			}
		})
	}

	go func() {
		wg.Wait()
		close(resultChan)
	}()

	select {
	case <-ctx.Done():
		return nil
	case result, ok := <-resultChan:
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
}

// mulAndSum computes result = matrix · vector in place for an integer
// coefficient vector (used to score candidate combinations quickly).
func mulAndSum(matrix *mat.Dense, vector []int, result []float64, rows int, cols int) {
	for row := range rows {
		result[row] = 0
		for col := range cols {
			result[row] += matrix.At(row, col) * float64(vector[col])
		}
	}
}

// computePseudoinverse returns the Moore-Penrose pseudoinverse of matrix as
// a cols×rows dense matrix. It is computed from the truncated SVD of the
// matrix: only the rank singular values above tol are kept (rank is the
// number of singular values greater than tol), which is the least-squares
// solution gonum's SolveTo expects. It fails if the factorization fails or
// the matrix has rank 0.
func computePseudoinverse(matrix *mat.Dense, tol float64) (*mat.Dense, error) {
	rows, cols := matrix.Dims()

	rank, svd, err := matrixRank(matrix, tol)
	if err != nil {
		return nil, err
	}

	if rank < 1 {
		return nil, fmt.Errorf("rank %d out of range", rank)
	}

	b := mat.NewDense(rows, rows, nil)
	for i := range rows {
		b.Set(i, i, 1.0)
	}

	inverse := mat.NewDense(cols, rows, nil)

	svd.SolveTo(inverse, b, rank)

	return inverse, nil
}

// matrixRank returns the rank of m — the number of singular values greater
// than tol — together with the SVD factorization, so callers that also need
// the decomposition can reuse it instead of factoring twice. An error is
// returned when SVD factorization fails.
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

// findNonZeroRows returns the indices of the rows of m that contain at least
// one element with an absolute value greater than tol. It drops the all-zero
// rows the [balancingAlgos.invAlgorithm] augmentation may produce.
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
