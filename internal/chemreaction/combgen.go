package chemreaction

import (
	"context"
	"fmt"
	"runtime"
	"sync"
)

type multiCombinationGenerator struct {
	maxCoef int
	k       int
}

func newMultiCombinationGenerator(maxCoef, k int) *multiCombinationGenerator {
	return &multiCombinationGenerator{
		maxCoef: maxCoef,
		k:       k,
	}
}

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
