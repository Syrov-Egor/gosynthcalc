package chemreaction

import (
	"context"
	"fmt"
	"runtime"
	"sync"
)

// multiCombinationGenerator enumerates candidate coefficient vectors: all
// combinations of k coefficients whose values range from 1 to maxCoef, in
// order of increasing maximal entry so that small solutions are found first.
type multiCombinationGenerator struct {
	maxCoef int // largest coefficient value to enumerate
	k       int // number of coefficients (compounds) in each combination
}

// newMultiCombinationGenerator returns a generator of combinations of k
// coefficients bounded by maxCoef.
func newMultiCombinationGenerator(maxCoef, k int) *multiCombinationGenerator {
	return &multiCombinationGenerator{
		maxCoef: maxCoef,
		k:       k,
	}
}

// generate streams every combination to the returned channel, spreading the
// work over numWorkers goroutines (runtime.GOMAXPROCS(0) when numWorkers is
// not positive). The channel is closed when the enumeration completes or ctx
// is cancelled. Progress ("Processing coef i of maxCoef") is printed to
// stdout while the search runs.
func (m *multiCombinationGenerator) generate(ctx context.Context, numWorkers int) <-chan []int {
	if numWorkers <= 0 {
		numWorkers = runtime.GOMAXPROCS(0)
	}

	out := make(chan []int, numWorkers*32)
	sem := make(chan struct{}, numWorkers)
	var wg sync.WaitGroup

	go func() {
		defer close(out)

		for maxValue := 1; maxValue <= m.maxCoef; maxValue++ {
			// Check for cancellation
			select {
			case <-ctx.Done():
				return
			default:
			}

			fmt.Printf("\r\033[2KProcessing coef %d of %d", maxValue, m.maxCoef)
			taskCount := maxValue
			workQueue := make(chan int, taskCount)

			for i := 1; i <= maxValue; i++ {
				select {
				case workQueue <- i:
				case <-ctx.Done():
					close(workQueue)
					return
				}
			}
			close(workQueue)

			var maxValueWg sync.WaitGroup
			maxValueWg.Add(taskCount)

			for i := 0; i < numWorkers; i++ {
				wg.Add(1)
				go func() {
					defer wg.Done()

					for startingValue := range workQueue {
						// Check for cancellation
						select {
						case <-ctx.Done():
							maxValueWg.Done()
							return
						default:
						}

						select {
						case sem <- struct{}{}:
						case <-ctx.Done():
							maxValueWg.Done()
							return
						}

						generateCombinations(ctx, startingValue, maxValue, m.k, out)
						<-sem
						maxValueWg.Done()
					}
				}()
			}
			maxValueWg.Wait()
		}
		wg.Wait()
	}()

	return out
}

// generateCombinations sends to out every k-length combination of values in
// [1, maxVal] whose first entry runs from startingValue up to maxVal while
// the remaining entries count like an odometer (the last one fastest); the
// enumeration therefore starts at (startingValue, 1, ..., 1). Each
// combination is copied before it is sent, so consumers may keep it. The
// function returns early when ctx is cancelled.
//
// Splitting the range of the first entry between callers is how
// [multiCombinationGenerator.generate] partitions the search space between
// its workers.
func generateCombinations(ctx context.Context, startingValue, maxVal, k int, out chan<- []int) {
	current := make([]int, k)
	current[0] = startingValue

	for i := 1; i < k; i++ {
		current[i] = 1
	}

	for {
		// Check for cancellation before processing
		select {
		case <-ctx.Done():
			return
		default:
		}

		result := make([]int, k)
		copy(result, current)

		// Send combination with cancellation support
		select {
		case out <- result:
		case <-ctx.Done():
			return
		}

		j := k - 1
		for j >= 0 && current[j] == maxVal {
			j--
		}

		if j < 0 || (j == 0 && current[0] > startingValue) {
			break
		}

		current[j]++

		for i := j + 1; i < k; i++ {
			current[i] = 1
		}
	}
}
