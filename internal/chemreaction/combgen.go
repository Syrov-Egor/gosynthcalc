package chemreaction

import (
	"fmt"
	"runtime"
	"sync"
)

type MultiCombinationGenerator struct {
	maxCoef int
	k       int
}

func NewMultiCombinationGenerator(maxCoef, k int) *MultiCombinationGenerator {
	return &MultiCombinationGenerator{
		maxCoef: maxCoef,
		k:       k,
	}
}

func (m *MultiCombinationGenerator) Generate(numWorkers int) <-chan []int {
	if numWorkers <= 0 {
		numWorkers = runtime.GOMAXPROCS(0)
	}

	out := make(chan []int, numWorkers*32)

	sem := make(chan struct{}, numWorkers)
	var wg sync.WaitGroup

	go func() {
		defer close(out)

		for maxValue := 1; maxValue <= m.maxCoef; maxValue++ {
			fmt.Printf("\r\033[2KProcessing coef %d of %d", maxValue, m.maxCoef)
			taskCount := maxValue
			workQueue := make(chan int, taskCount)

			for i := 1; i <= maxValue; i++ {
				workQueue <- i
			}
			close(workQueue)
			var maxValueWg sync.WaitGroup
			maxValueWg.Add(taskCount)

			for range numWorkers {
				wg.Add(1)
				go func() {
					defer wg.Done()

					for startingValue := range workQueue {
						sem <- struct{}{}
						generateCombinations(startingValue, maxValue, m.k, out)
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

func generateCombinations(startingValue, maxVal, k int, out chan<- []int) {
	current := make([]int, k)
	current[0] = startingValue

	for i := 1; i < k; i++ {
		current[i] = 1
	}

	for {
		result := make([]int, k)
		copy(result, current)
		out <- result
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
