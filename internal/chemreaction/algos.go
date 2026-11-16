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

// nullSpaceByInversion computes a basis of the null space of the reaction
// matrix using Thorne's matrix-inversion trick
// (https://arxiv.org/abs/1110.4321).
//
// Let A be the m×n composition matrix and r its rank. The idea is to pick r
// linearly independent rows R of A (r×n), to append n-r unit rows e_f (one
// per "free" column f) so that the stacked matrix
//
//	M = [ R ]
//	    [ E ]
//
// is square and invertible, and to read the null space off the last n-r
// columns of M⁻¹. Indeed M·M⁻¹ = I, so for the k-th of those columns x_k:
// R·x_k = 0 (hence A·x_k = 0, the other rows of A being combinations of the
// rows of R) and e_{f_j}·x_k = δ_jk.
//
// Steps:
//
//  1. The rank r and the nullity n-r are supplied by the balancer, which
//     determines them once from the singular values of A.
//  2. A Gaussian elimination with partial pivoting, left to right, tracks
//     which ORIGINAL rows end up as pivot rows (they are independent, which
//     drops dependent rows such as a charge row that is a combination of the
//     element rows) and which columns carry a pivot. The columns without a
//     pivot are the free columns. For a generic matrix these are the
//     rightmost columns, exactly as in the article, but unlike a fixed
//     "identity in the last columns" the choice always yields an invertible M.
//  3. M is built from the original rows and the unit rows, inverted, and its
//     last n-r columns are returned (an n×(n-r) matrix; column k is the k-th
//     basis vector, indexed by compound).
//
// The original rows, and not an orthonormal basis of the row space, are used
// on purpose: M then stays as exact as the input. For H + O₂ = H₁₀₀₀₀₀₁O₂ the
// matrix [[1 0 N] [0 2 2]] has a condition number of about N/2, so an SVD
// basis (the first r rows of Vᵀ) loses about log₁₀(cond) digits and the
// coefficient N comes out as 1000001.0003, while with the integer rows the
// augmented matrix is [[1 0 N] [0 2 2] [0 0 1]] and its inverse is exact.
//
// Safeguards: every returned basis vector is verified (scaled residual of
// A·x). If the exact-rows route fails or is inaccurate, which can only happen
// for badly conditioned matrices that no chemical equation produces, the same
// construction is repeated with the orthonormal SVD basis, which is
// backward-stable but not exact. If that fails the check as well, an error is
// returned: a wrong null space is never returned silently.
//
// The rank and nullity of the reaction matrix are supplied by the caller
// (see [balancer.reactionRankAndNullity]) so that they are not computed more
// than once. The basis is nil when the nullity is zero.
func (b *balancingAlgos) nullSpaceByInversion(rank, nullity int) (basis *mat.Dense, err error) {
	m, n := b.ReactionMatrix.Dims()
	if n == 0 || m == 0 {
		return nil, fmt.Errorf("empty reaction matrix")
	}

	if nullity == 0 {
		return nil, nil
	}
	if rank == 0 {
		// The zero matrix: every unit vector is a null vector.
		basis = mat.NewDense(n, n, nil)
		for i := range n {
			basis.Set(i, i, 1)
		}
		return basis, nil
	}

	// 1. Exact route: independent ORIGINAL rows.
	var exactErr error
	pivotRows, free, exactErr := pivotRowsAndFreeColumns(b.ReactionMatrix, rank, b.Tolerance)
	if exactErr == nil && len(free) == nullity {
		rows := mat.NewDense(rank, n, nil)
		for i, r := range pivotRows {
			for j := range n {
				rows.Set(i, j, b.ReactionMatrix.At(r, j))
			}
		}
		basis, exactErr = augmentedInverseBasis(rows, free)
		if exactErr == nil && !b.nullBasisIsAccurate(basis) {
			exactErr = fmt.Errorf("inaccurate null space")
		}
		if exactErr == nil {
			return basis, nil
		}
	} else if exactErr == nil {
		exactErr = fmt.Errorf("rank/pivot mismatch: expected %d free columns, found %d", nullity, len(free))
	}

	// 2. Fallback: orthonormal basis of the row space (first rank rows of Vᵀ).
	var svd mat.SVD
	if !svd.Factorize(b.ReactionMatrix, mat.SVDThin) {
		return nil, fmt.Errorf("SVD factorization failed: %w", exactErr)
	}
	var v mat.Dense
	svd.VTo(&v)
	rowSpace := mat.NewDense(rank, n, nil)
	for i := range rank {
		for j := range n {
			rowSpace.Set(i, j, v.At(j, i))
		}
	}
	_, free, err = pivotRowsAndFreeColumns(rowSpace, rank, b.Tolerance)
	if err == nil && len(free) == nullity {
		basis, err = augmentedInverseBasis(rowSpace, free)
		if err == nil && b.nullBasisIsAccurate(basis) {
			return basis, nil
		}
	}

	return nil, fmt.Errorf("could not compute an accurate null space (the reaction matrix is numerically ill-conditioned): %w", exactErr)
}

// augmentedInverseBasis stacks rows (rank×n) on top of the unit rows of the
// free columns, inverts the resulting n×n matrix and returns its last
// n-rank columns, which form a basis of the null space of rows.
func augmentedInverseBasis(rows *mat.Dense, free []int) (*mat.Dense, error) {
	rank, n := rows.Dims()

	mm := mat.NewDense(n, n, nil)
	mm.Slice(0, rank, 0, n).(*mat.Dense).Copy(rows)
	for k, f := range free {
		mm.Set(rank+k, f, 1)
	}

	var inv mat.Dense
	if err := inv.Inverse(mm); err != nil {
		return nil, fmt.Errorf("matrix inversion failed: %w", err)
	}

	nullity := n - rank
	basis := mat.NewDense(n, nullity, nil)
	for k := range nullity {
		for i := range n {
			basis.Set(i, k, inv.At(i, rank+k))
		}
	}
	return basis, nil
}

// nullBasisIsAccurate checks that every column x of basis satisfies A·x ≈ 0.
// For each row i the residual |Σ a_ij·x_j| must be tiny compared to the terms
// that were added up (Σ |a_ij·x_j|, a component-wise test), or, to tolerate
// the rounding noise that sits in entries that are really zero, tiny compared
// to the size of the row times the size of x (a norm-wise test).
func (b *balancingAlgos) nullBasisIsAccurate(basis *mat.Dense) bool {
	m, n := b.ReactionMatrix.Dims()
	_, k := basis.Dims()
	for c := range k {
		xMax := 0.0
		for j := range n {
			xMax = math.Max(xMax, math.Abs(basis.At(j, c)))
		}
		for i := range m {
			var sum, scale, rowAbs float64
			for j := range n {
				t := b.ReactionMatrix.At(i, j) * basis.At(j, c)
				sum += t
				scale += math.Abs(t)
				rowAbs += math.Abs(b.ReactionMatrix.At(i, j))
			}
			if math.Abs(sum) > 1e-10*scale+1e-13*rowAbs*xMax {
				return false
			}
		}
	}
	return true
}

// pivotRowsAndFreeColumns runs a left-to-right Gaussian elimination with
// partial pivoting on a copy of a and returns
//
//   - the indices (ascending) of the ORIGINAL rows of a that were used as
//     pivot rows: they are linearly independent and span the row space;
//   - the indices of the columns without a pivot (the free columns).
//
// The stacked matrix [those rows; unit rows of the free columns] is
// invertible: the pivot rows are, up to a unit lower-triangular transform,
// the echelon rows, whose pivot-column minor is triangular and nonsingular.
// It returns an error if the number of pivots differs from rank.
func pivotRowsAndFreeColumns(a *mat.Dense, rank int, tol float64) (pivotRows, free []int, err error) {
	m, n := a.Dims()
	w := mat.DenseCopyOf(a)

	maxAbs := 0.0
	for i := range m {
		for j := range n {
			maxAbs = math.Max(maxAbs, math.Abs(w.At(i, j)))
		}
	}
	eps := tol * math.Max(1, maxAbs)

	orig := make([]int, m) // orig[i] = original index of the row now at i
	for i := range orig {
		orig[i] = i
	}
	isPivot := make([]bool, n)

	row := 0
	for col := 0; col < n && row < m; col++ {
		best, bestVal := -1, eps
		for i := row; i < m; i++ {
			if v := math.Abs(w.At(i, col)); v > bestVal {
				best, bestVal = i, v
			}
		}
		if best < 0 {
			continue
		}
		if best != row {
			for j := range n {
				t := w.At(row, j)
				w.Set(row, j, w.At(best, j))
				w.Set(best, j, t)
			}
			orig[row], orig[best] = orig[best], orig[row]
		}
		piv := w.At(row, col)
		for i := row + 1; i < m; i++ {
			f := w.At(i, col) / piv
			if f == 0 {
				continue
			}
			for j := col; j < n; j++ {
				w.Set(i, j, w.At(i, j)-f*w.At(row, j))
			}
		}
		isPivot[col] = true
		row++
	}
	if row != rank {
		return nil, nil, fmt.Errorf("could not find %d pivot rows (found %d)", rank, row)
	}

	pivotRows = slices.Clone(orig[:row])
	slices.Sort(pivotRows)
	for j, p := range isPivot {
		if !p {
			free = append(free, j)
		}
	}
	return pivotRows, free, nil
}

// invAlgorithm computes the coefficients with the reaction-matrix
// inverse algorithm proposed by Thorne.
//
// It is a thin layer over [balancingAlgos.nullSpaceByInversion]:
//
//   - nullity 0: no balanced equation exists → error;
//   - nullity ≥ 2: there are infinitely many independent balanced equations,
//     so there is no single answer → error (the basis is still available
//     from nullSpaceByInversion);
//   - nullity 1: the null vector is unique up to scale. It is oriented so
//     that the reactants are positive, validated (every compound must take
//     part, and reactants and products must have opposite signs), and then
//     divided by its smallest magnitude. Rounding to integers is left to the
//     caller, as before.
//
// The reaction matrix is assumed to hold plain (non-negative) atom counts for
// both sides, as in the article, so a balanced vector has opposite signs on
// the two sides.
//
// rank and nullity are supplied by the balancer; they are not recomputed
// here.
func (b *balancingAlgos) invAlgorithm(rank, nullity int) ([]float64, error) {
	_, cols := b.ReactionMatrix.Dims()
	if b.SeparatorPos <= 0 || b.SeparatorPos >= cols {
		return nil, fmt.Errorf("reaction needs at least one reactant and one product")
	}

	switch {
	case nullity == 0:
		return nil, fmt.Errorf("the reaction cannot be balanced (nullity 0)")
	case nullity > 1:
		return nil, fmt.Errorf("the reaction has nullity %d: it has infinitely many independent balanced equations, not a single one", nullity)
	}

	basis, err := b.nullSpaceByInversion(rank, nullity)
	if err != nil {
		return nil, err
	}

	x := mat.Col(nil, 0, basis)

	maxAbs := 0.0
	for _, val := range x {
		maxAbs = math.Max(maxAbs, math.Abs(val))
	}
	eps := b.Tolerance * maxAbs

	// Orient the vector so that reactants come out positive.
	if x[0] < 0 {
		for i := range x {
			x[i] = -x[i]
		}
	}

	for i, val := range x {
		switch {
		case math.Abs(val) <= eps:
			return nil, fmt.Errorf("compound #%d does not take part in the only balanced equation", i)
		case i < b.SeparatorPos && val < 0:
			return nil, fmt.Errorf("compound #%d (reactant) has the wrong sign: it belongs on the product side", i)
		case i >= b.SeparatorPos && val > 0:
			return nil, fmt.Errorf("compound #%d (product) has the wrong sign: it belongs on the reactant side", i)
		}
	}

	coefs := make([]float64, len(x))
	for i, val := range x {
		coefs[i] = math.Abs(val)
	}
	minVal := slices.Min(coefs)
	for i := range coefs {
		coefs[i] /= minVal
	}
	return coefs, nil
}

// gPInvAlgorithm computes the coefficients with the general pseudoinverse
// algorithm proposed by
// [Risteski](http://koreascience.or.kr/article/JAKO201314358624990.page).
// The algorithm can be described in steps:
//
//  1. Stack the reactant matrix and the negative product matrix: the
//     balanced coefficients are exactly the solutions of A·x = 0.
//
//  2. Calculate the MP pseudoinverse A⁺ of this matrix.
//
//  3. Calculate the coefficients by the general solution of a homogeneous
//     system: x = (I – A⁺A)a, where I is the identity matrix and a is an
//     arbitrary vector (here a vector of ones).
//
// (I – A⁺A) is the orthogonal projector onto the null space of A, so the
// result is the projection of the vector of ones onto the set of balanced
// solutions. The consequences:
//
//   - nullity 0: the only solution is x = 0, the reaction cannot be balanced;
//     this is reported as an error instead of returning float noise;
//   - nullity 1: the projection is a positive multiple of the unique balanced
//     vector, so the result is already positive if the reaction is feasible as
//     written (mixed signs mean that some compound is on the wrong side);
//   - nullity ≥ 2: the result is just one of infinitely many solutions, the
//     one closest to a vector of ones. Nothing guarantees that it is positive
//     or free of zeros even when a fully positive solution exists, so the
//     caller must validate it (and fall back to another search if needed).
//
// The result does not depend on the order of the rows (elements) of the
// reaction matrix: A⁺A is the orthogonal projector onto the row space, and
// permuting rows does not change that space. Only floating point noise
// changes.
//
// The output is a float array, so to get integer coefficients it has to be
// converted, which is done (together with validation) by
// [balancer.calculateByMethod].
//
// rank and nullity are supplied by the balancer; they are not recomputed
// here.
func (b *balancingAlgos) gPInvAlgorithm(rank, nullity int) ([]float64, error) {
	rows, cols := b.ReactionMatrix.Dims()

	matrix := mat.NewDense(rows, cols, nil)
	matrix.Copy(b.ReactionMatrix)

	for i := range rows {
		for j := b.SeparatorPos; j < cols; j++ {
			matrix.Set(i, j, -matrix.At(i, j))
		}
	}

	if rank == 0 {
		return nil, fmt.Errorf("the reaction matrix contains no atoms")
	}
	if nullity == 0 {
		return nil, fmt.Errorf("the reaction cannot be balanced (nullity 0)")
	}

	inverse, err := pseudoinverseWithRank(matrix, rank)
	if err != nil {
		return nil, err
	}

	// x = (I - A⁺A)·1
	temp := mat.NewDense(cols, cols, nil)
	temp.Product(inverse, matrix)
	temp.Scale(-1, temp)
	for i := range cols {
		temp.Set(i, i, temp.At(i, i)+1)
	}

	ones := mat.NewVecDense(cols, nil)
	for i := range cols {
		ones.SetVec(i, 1)
	}
	coefs := mat.NewVecDense(cols, nil)
	coefs.MulVec(temp, ones)

	result := make([]float64, cols)
	for i := range cols {
		result[i] = coefs.AtVec(i)
	}
	zeroNoise(result, b.Tolerance)

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
// G = 0 is a normal case (it happens whenever every product column lies in
// the column space of A, e.g. when A has full row rank, as in
// H₂ + O₂ = H₂O), and then G⁺ = 0 and y = u. The pseudoinverse of a zero
// matrix is therefore the zero matrix and not an error.
//
// The nullity of the whole reaction matrix is
// (reactants − rank A) + (products − rank G), so it comes for free; nullity 0
// is reported as an error. As with [balancingAlgos.gPInvAlgorithm], for
// nullity ≥ 2 the result is only one particular solution and the caller has
// to validate it. Note that the term (I-A^-A)v is an arbitrary element of the
// null space of A (reactants that are dependent on each other), so when that
// space is not trivial the result can have wrong signs.
//
// rank, nullity and reactantRank are supplied by the balancer. Since
// rank G = rank − rank A, the rank of the intermediate G matrix is derived
// from them instead of being computed again with an SVD.
//
// Note: while this algorithm and [balancingAlgos.gPInvAlgorithm] are very
// similar, there are some differences in output results. This method exists
// mostly for legacy purposes, like balancing some reactions according to
// Risteski.
func (b *balancingAlgos) pPInvAlgorithm(rank, nullity, reactantRank int) ([]float64, error) {
	reactantRows, reactantCols := b.ReactantMatrix.Dims()

	if nullity == 0 {
		return nil, fmt.Errorf("the reaction cannot be balanced (nullity 0)")
	}

	mpInverse, err := pseudoinverseWithRank(b.ReactantMatrix, reactantRank)
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

	gRank := rank - reactantRank
	if gRank < 0 || gRank > min(reactantRows, b.ProductCols) {
		return nil, fmt.Errorf("invalid rank of the G matrix: %d", gRank)
	}
	gPinv, err := pseudoinverseWithRank(mat.DenseCopyOf(&gMatrix), gRank)
	if err != nil {
		return nil, fmt.Errorf("error computing pseudoinverse of the G matrix: %w", err)
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
	zeroNoise(coefs, b.Tolerance)

	return coefs, nil
}

// pseudoinverseWithRank returns the Moore-Penrose pseudoinverse of matrix as
// a cols×rows dense matrix using the supplied rank. The rank is determined by
// the balancer and is not recomputed here.
//
// It is computed from the truncated SVD of the matrix: only the first rank
// singular values are kept, which is the least-squares solution gonum's
// SolveTo expects. A zero rank (in particular for a zero matrix) gives the
// zero matrix as its pseudoinverse; gonum's SolveTo cannot handle rank 0, so
// that case is returned directly. The function fails only if the
// factorization fails.
func pseudoinverseWithRank(matrix *mat.Dense, rank int) (*mat.Dense, error) {
	rows, cols := matrix.Dims()

	inverse := mat.NewDense(cols, rows, nil)
	if rank < 1 {
		return inverse, nil
	}

	var svd mat.SVD
	if !svd.Factorize(matrix, mat.SVDFull) {
		return nil, fmt.Errorf("SVD factorization failed")
	}

	b := mat.NewDense(rows, rows, nil)
	for i := range rows {
		b.Set(i, i, 1.0)
	}

	svd.SolveTo(inverse, b, rank)

	return inverse, nil
}

// zeroNoise sets to zero the entries of v that are negligible compared to the
// largest one (floating point residue of a projection onto the null space),
// so that "this compound does not take part" is a clean 0 for the caller.
func zeroNoise(v []float64, tol float64) {
	maxAbs := 0.0
	for _, x := range v {
		maxAbs = math.Max(maxAbs, math.Abs(x))
	}
	eps := tol * math.Max(1, maxAbs)
	for i, x := range v {
		if math.Abs(x) <= eps {
			v[i] = 0
		}
	}
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
