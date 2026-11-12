package chemreaction

import (
	"context"
	"runtime"
	"sync"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

// combinatorial finds a solution of a Diophantine matrix equation by
// enumerating candidate coefficient vectors in [1,maxCoef]^n. It returns the
// lexicographically first valid vector within the bound, so the result is
// deterministic regardless of scheduling.
//
// The lexicographic search space is partitioned into contiguous ranges of the
// first coefficient. Those ranges are searched concurrently by a bounded
// worker pool, and results are consumed in range order: a later range can run
// speculatively, but it is accepted only after every earlier range has proved
// to contain no solution. This preserves exactly the sequential lexicographic
// answer while using multiple cores. The search is silent and honours ctx.
//
// The method is incomplete because the true solution may contain a coefficient
// above maxCoef, and it is exponential in the number of compounds. It is a
// bounded diagnostic and is not used by [balancer.Auto].
func (b *balancingAlgos) combinatorial(ctx context.Context, maxCoef uint) []float64 {
	return b.combinatorialSearch(ctx, maxCoef, true)
}

// combinatorialSequential is the single-goroutine reference implementation of
// [balancingAlgos.combinatorial]. It returns the same vector and is used by
// tests and benchmarks to verify the parallel search.
func (b *balancingAlgos) combinatorialSequential(ctx context.Context, maxCoef uint) []float64 {
	return b.combinatorialSearch(ctx, maxCoef, false)
}

func (b *balancingAlgos) combinatorialSearch(ctx context.Context, maxCoef uint, parallel bool) []float64 {
	if ctx == nil {
		ctx = context.Background()
	}
	if ctx.Err() != nil {
		return nil
	}

	_, cols := b.ReactionMatrix.Dims()
	if maxCoef == 0 || cols == 0 {
		return nil
	}
	limit := int(maxCoef)
	current := make([]int, cols)

	var (
		found []int
		ok    bool
	)
	if parallel {
		found, ok = b.searchLexParallel(ctx, current, limit)
	} else {
		found, ok = b.searchLexSequential(ctx, current, 0, limit)
	}
	if !ok {
		return nil
	}

	result := make([]float64, len(found))
	for i, v := range found {
		result[i] = float64(v)
	}
	return result
}

// searchLexSequential fills current[pos:] with the lexicographically first
// solution in [1,limit]^(len(current)-pos). Prefix entries before pos are
// fixed by the caller.
func (b *balancingAlgos) searchLexSequential(ctx context.Context, current []int, pos, limit int) ([]int, bool) {
	if ctx.Err() != nil {
		return nil, false
	}

	reacSum := make([]float64, b.ReactantRows)
	prodSum := make([]float64, b.ProductRows)

	var visit func(int) ([]int, bool)
	visit = func(p int) ([]int, bool) {
		if ctx.Err() != nil {
			return nil, false
		}
		if p == len(current) {
			if b.balancesCandidate(current, reacSum, prodSum) {
				return append([]int(nil), current...), true
			}
			return nil, false
		}
		for v := 1; v <= limit; v++ {
			current[p] = v
			if solution, ok := visit(p + 1); ok {
				return solution, true
			}
			if ctx.Err() != nil {
				return nil, false
			}
		}
		return nil, false
	}
	return visit(pos)
}

// searchLexRangeSequential tries values lo..hi at current[pos] in order and
// otherwise performs the same lexicographic search as
// [balancingAlgos.searchLexSequential].
func (b *balancingAlgos) searchLexRangeSequential(ctx context.Context, current []int, pos, lo, hi, limit int) ([]int, bool) {
	for v := lo; v <= hi; v++ {
		if ctx.Err() != nil {
			return nil, false
		}
		current[pos] = v
		if solution, ok := b.searchLexSequential(ctx, current, pos+1, limit); ok {
			return solution, true
		}
	}
	return nil, false
}

// combParallelWorkers returns the number of concurrently active branch
// searches. Long-running CPU-bound branches scale best with a modest number of
// workers: oversubscribing simultaneous multithreading threads makes the
// critical lexicographic branch slower. Small machines still get at least
// four workers (or all available pseudo-processors when there are fewer).
func combParallelWorkers() int {
	procs := runtime.GOMAXPROCS(0)
	if procs < 1 {
		return 1
	}
	workers := procs / 4
	if workers < 4 {
		workers = 4
	}
	if workers > procs {
		workers = procs
	}
	return workers
}

// searchLexParallel searches the same lexicographic space as
// [balancingAlgos.searchLexSequential], but partitions the first coefficient
// into contiguous ranges and searches those ranges with a bounded worker pool.
// Results are consumed in range order, so the returned vector is the same as
// the sequential search.
func (b *balancingAlgos) searchLexParallel(ctx context.Context, current []int, limit int) ([]int, bool) {
	workers := combParallelWorkers()
	if workers <= 1 {
		return b.searchLexSequential(ctx, current, 0, limit)
	}

	// Keep the branch count bounded while still separating the first
	// coefficient values for every realistic maxCoef bound. Each branch is a
	// contiguous lexicographic range; the worker count limits how many of
	// them run at once.
	const maxCombBranches = 256
	branchCount := min(limit, maxCombBranches)
	if branchCount < 1 {
		branchCount = 1
	}

	type branchResult struct {
		solution []int
		ok       bool
	}

	childCtx, cancel := context.WithCancel(ctx)
	defer cancel()

	sem := make(chan struct{}, workers)
	branches := make([]chan branchResult, branchCount)
	var wg sync.WaitGroup

	launch := func(i int, slotHeld bool) {
		ch := make(chan branchResult, 1)
		branches[i] = ch
		wg.Add(1)
		go func() {
			defer wg.Done()
			if slotHeld {
				defer func() { <-sem }()
			} else {
				select {
				case sem <- struct{}{}:
					defer func() { <-sem }()
				case <-childCtx.Done():
					ch <- branchResult{}
					return
				}
			}

			lo := i*limit/branchCount + 1
			hi := (i + 1) * limit / branchCount
			branchCurrent := make([]int, len(current))
			copy(branchCurrent, current)

			solution, ok := b.searchLexRangeSequential(childCtx, branchCurrent, 0, lo, hi, limit)
			ch <- branchResult{solution: solution, ok: ok}
		}()
	}

	// Start the earliest ranges first and reserve their worker slots so the
	// lexicographically critical work cannot be starved by speculative later
	// ranges.
	initial := min(workers, branchCount)
	for i := 0; i < initial; i++ {
		sem <- struct{}{}
		launch(i, true)
	}
	for i := initial; i < branchCount; i++ {
		launch(i, false)
	}

	for i := 0; i < branchCount; i++ {
		result := <-branches[i]
		if result.ok {
			cancel()
			wg.Wait()
			return result.solution, true
		}
	}

	wg.Wait()
	return nil, false
}

// balancesCandidate reports whether the integer candidate vector balances the
// reaction within Tolerance. reacSum and prodSum are scratch buffers reused by
// the sequential search.
func (b *balancingAlgos) balancesCandidate(candidate []int, reacSum, prodSum []float64) bool {
	reactantCoefs := candidate[:b.SeparatorPos]
	productCoefs := candidate[b.SeparatorPos:]
	mulAndSum(b.ReactantMatrix, reactantCoefs, reacSum, b.ReactantRows, b.ReactantCols)
	mulAndSum(b.ProductMatrix, productCoefs, prodSum, b.ProductRows, b.ProductCols)
	return floats.EqualApprox(reacSum, prodSum, b.Tolerance)
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
