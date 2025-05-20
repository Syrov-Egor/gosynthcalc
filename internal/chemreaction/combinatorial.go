package chemreaction

import (
	"fmt"
	"sync"
)

// MultiCombinationGenerator produces all possible k-combinations with repetition
// where each element is between 1 and maxCoef with increasing maxCoef values
type MultiCombinationGenerator struct {
	maxCoef int
	k       int
}

// NewMultiCombinationGenerator creates a new generator
func NewMultiCombinationGenerator(maxCoef, k int) *MultiCombinationGenerator {
	return &MultiCombinationGenerator{
		maxCoef: maxCoef,
		k:       k,
	}
}

// GenerateForMaxValue produces combinations for a specific max value and k length
func (m *MultiCombinationGenerator) GenerateForMaxValue(maxValue, k int, workers int, out chan<- []int) {
	var wg sync.WaitGroup

	// Function to generate combinations with specific maxValue
	generateWorker := func(startingValue, maxVal, k int, wg *sync.WaitGroup) {
		defer wg.Done()

		// Initialize current combination
		current := make([]int, k)
		current[0] = startingValue

		// Setup initial state
		for i := 1; i < k; i++ {
			current[i] = 1
		}

		for {
			// Make a copy of the current combination and send it
			result := make([]int, k)
			copy(result, current)
			out <- result

			// Generate next combination in lexicographic order
			j := k - 1
			for j >= 0 && current[j] == maxVal {
				j--
			}

			// If we've exhausted all combinations starting with startingValue
			if j < 0 || (j == 0 && current[0] > startingValue) {
				break
			}

			current[j]++

			// Reset subsequent positions
			for i := j + 1; i < k; i++ {
				current[i] = 1
			}
		}
	}

	// Start workers based on the first position values
	for i := 1; i <= maxValue; i += workers {
		for w := 0; w < workers && i+w <= maxValue; w++ {
			wg.Add(1)
			go generateWorker(i+w, maxValue, k, &wg)
		}
		wg.Wait()
	}
}

// Generate produces all combinations with increasing max values from 1 to maxCoef
func (m *MultiCombinationGenerator) Generate(workers int) <-chan []int {
	out := make(chan []int)

	go func() {
		// For each maxValue from 1 to maxCoef
		for maxValue := 1; maxValue <= m.maxCoef; maxValue++ {

			m.GenerateForMaxValue(maxValue, m.k, workers, out)
			fmt.Printf("\r\033[2KCoef %d of %d", maxValue+1, m.maxCoef)
		}
		close(out)

	}()

	return out
}
