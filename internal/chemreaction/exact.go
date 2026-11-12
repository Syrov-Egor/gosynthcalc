package chemreaction

import (
	"fmt"
	"math/big"

	"gonum.org/v1/gonum/mat"
)

// exactMatrix is the signed reaction matrix A = [R | -P] over the rationals.
// Reactant columns are positive, product columns are negative.
type exactMatrix struct {
	rows, cols int
	data       [][]*big.Rat // data[i][j]
}

// matrixAnalysis is the cached result of exact row reduction of an
// [exactMatrix].
type matrixAnalysis struct {
	rank    int
	nullity int
	pivots  []int
	rref    [][]*big.Rat // rank rows, cols columns
	basis   [][]*big.Rat // cols rows, nullity columns
}

// newExactMatrix converts an unsigned float reaction matrix into the signed
// exact matrix A = [R | -P]. The first separatorPos columns are reactants and
// stay positive; the remaining product columns are negated.
//
// Each float entry is converted with [ratFromFloat]. NaN and infinities are
// rejected with [ErrNonFiniteMatrix]; a nil matrix, a non-positive column
// count or an out-of-range separator position is rejected with
// [ErrInvalidMatrix].
func newExactMatrix(m *mat.Dense, separatorPos int) (*exactMatrix, error) {
	if m == nil {
		return nil, fmt.Errorf("%w: nil reaction matrix", ErrInvalidMatrix)
	}
	rows, cols := m.Dims()
	if cols <= 0 {
		return nil, fmt.Errorf("%w: reaction matrix has %d columns", ErrInvalidMatrix, cols)
	}
	if separatorPos < 0 || separatorPos > cols {
		return nil, fmt.Errorf("%w: separator position %d is outside [0,%d]",
			ErrInvalidMatrix, separatorPos, cols)
	}

	data := make([][]*big.Rat, rows)
	for i := 0; i < rows; i++ {
		data[i] = make([]*big.Rat, cols)
		for j := 0; j < cols; j++ {
			v, err := ratFromFloat(m.At(i, j))
			if err != nil {
				return nil, fmt.Errorf("matrix entry (%d,%d): %w", i, j, err)
			}
			if j >= separatorPos {
				v.Neg(v)
			}
			data[i][j] = v
		}
	}
	return &exactMatrix{rows: rows, cols: cols, data: data}, nil
}

// validate checks that the matrix dimensions agree with its data slice.
func (A *exactMatrix) validate() error {
	if A == nil {
		return fmt.Errorf("%w: nil exact matrix", ErrInvalidMatrix)
	}
	if A.cols <= 0 {
		return fmt.Errorf("%w: exact matrix has %d columns", ErrInvalidMatrix, A.cols)
	}
	if A.rows < 0 || len(A.data) != A.rows {
		return fmt.Errorf("%w: exact matrix declares %d rows but stores %d",
			ErrInvalidMatrix, A.rows, len(A.data))
	}
	for i := range A.data {
		if len(A.data[i]) != A.cols {
			return fmt.Errorf("%w: exact matrix row %d has %d entries, want %d",
				ErrInvalidMatrix, i, len(A.data[i]), A.cols)
		}
		for j := range A.data[i] {
			if A.data[i][j] == nil {
				return fmt.Errorf("%w: exact matrix entry (%d,%d) is nil",
					ErrInvalidMatrix, i, j)
			}
		}
	}
	return nil
}

// exactRREF computes the reduced row-echelon form of A over the rationals.
// It returns the non-zero rows of the RREF, the pivot column of each row and
// the rank. The computation is deterministic: it always pivots on the leftmost
// non-zero column and the first non-zero row below the current one, and it
// uses exact [big.Rat] arithmetic with no tolerance.
func exactRREF(A *exactMatrix) ([][]*big.Rat, []int, int) {
	if A == nil || A.cols <= 0 || A.rows < 0 || len(A.data) != A.rows {
		return nil, nil, 0
	}

	m, n := A.rows, A.cols
	e := make([][]*big.Rat, m)
	for i := 0; i < m; i++ {
		if len(A.data[i]) != n {
			return nil, nil, 0
		}
		e[i] = make([]*big.Rat, n)
		for j := 0; j < n; j++ {
			e[i][j] = new(big.Rat).Set(A.data[i][j])
		}
	}

	row := 0
	pivots := make([]int, 0, min(m, n))
	for col := 0; col < n && row < m; col++ {
		p := -1
		for i := row; i < m; i++ {
			if e[i][col].Sign() != 0 {
				p = i
				break
			}
		}
		if p < 0 {
			continue
		}
		if p != row {
			e[row], e[p] = e[p], e[row]
		}

		pivot := new(big.Rat).Set(e[row][col])
		for j := col; j < n; j++ {
			e[row][j].Quo(e[row][j], pivot)
		}

		for i := 0; i < m; i++ {
			if i == row || e[i][col].Sign() == 0 {
				continue
			}
			f := new(big.Rat).Set(e[i][col])
			for j := col; j < n; j++ {
				e[i][j].Sub(e[i][j], new(big.Rat).Mul(f, e[row][j]))
			}
		}

		pivots = append(pivots, col)
		row++
	}

	return e[:row], pivots, row
}

// nullSpaceBasis builds a basis of the null space from an RREF and its pivot
// columns. If the free columns are f_0, ..., f_{k-1}, the j-th basis vector
// has 1 in coordinate f_j and -RREF[i][f_j] in pivot coordinate pivots[i].
// The result has one row per column of the original matrix and one column per
// basis vector (the nullity). A malformed RREF returns nil.
func nullSpaceBasis(rref [][]*big.Rat, pivots []int, n int) [][]*big.Rat {
	rank := len(rref)
	if n <= 0 || rank > n || len(pivots) != rank {
		return nil
	}
	k := n - rank
	if k <= 0 {
		return nil
	}

	isPivot := make([]bool, n)
	for _, c := range pivots {
		if c < 0 || c >= n || isPivot[c] {
			return nil
		}
		isPivot[c] = true
	}
	free := make([]int, 0, k)
	for c := 0; c < n; c++ {
		if !isPivot[c] {
			free = append(free, c)
		}
	}
	if len(free) != k {
		return nil
	}
	for i := range rref {
		if len(rref[i]) != n {
			return nil
		}
	}

	basis := make([][]*big.Rat, n)
	for i := 0; i < n; i++ {
		basis[i] = make([]*big.Rat, k)
		for j := 0; j < k; j++ {
			basis[i][j] = new(big.Rat)
		}
	}
	for j, f := range free {
		basis[f][j].SetInt64(1)
	}
	for i := 0; i < rank; i++ {
		for j, f := range free {
			basis[pivots[i]][j].Neg(rref[i][f])
		}
	}
	return basis
}

// exactDenseRank returns the exact rank of a float matrix, treating all
// entries as positive (no product-column negation). It is used by the legacy
// pseudoinverse heuristics to choose their SVD truncation rank from exact
// arithmetic instead of an absolute tolerance.
func exactDenseRank(m *mat.Dense) (int, error) {
	if m == nil {
		return 0, fmt.Errorf("%w: nil matrix", ErrInvalidMatrix)
	}
	_, cols := m.Dims()
	if cols <= 0 {
		return 0, fmt.Errorf("%w: matrix has %d columns", ErrInvalidMatrix, cols)
	}
	A, err := newExactMatrix(m, cols)
	if err != nil {
		return 0, err
	}
	a, err := analyze(A)
	if err != nil {
		return 0, err
	}
	return a.rank, nil
}

// analyze validates A, computes its exact RREF, rank, nullity and null-space
// basis, and returns them together. The null-space basis is nil when the
// nullity is zero; otherwise it has A.cols rows and A.nullity columns and
// satisfies A * basis = 0 exactly.
func analyze(A *exactMatrix) (matrixAnalysis, error) {
	if err := A.validate(); err != nil {
		return matrixAnalysis{}, err
	}

	rref, pivots, rank := exactRREF(A)
	if len(pivots) != rank {
		return matrixAnalysis{}, fmt.Errorf("%w: inconsistent row reduction", ErrInvalidMatrix)
	}
	nullity := A.cols - rank
	if nullity < 0 {
		return matrixAnalysis{}, fmt.Errorf("%w: rank %d exceeds column count %d",
			ErrInvalidMatrix, rank, A.cols)
	}

	var basis [][]*big.Rat
	if nullity > 0 {
		basis = nullSpaceBasis(rref, pivots, A.cols)
		if basis == nil {
			return matrixAnalysis{}, fmt.Errorf("%w: failed to build null-space basis", ErrInvalidMatrix)
		}
	}

	return matrixAnalysis{
		rank:    rank,
		nullity: nullity,
		pivots:  pivots,
		rref:    rref,
		basis:   basis,
	}, nil
}
