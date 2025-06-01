package utils

import (
	"math"
	"math/big"
	"reflect"
	"testing"
)

func TestStringCounter(t *testing.T) {
	tests := []struct {
		name     string
		input    string
		expected map[string]int
	}{
		{
			name:     "simple string",
			input:    "hello",
			expected: map[string]int{"h": 1, "e": 1, "l": 2, "o": 1},
		},
		{
			name:     "empty string",
			input:    "",
			expected: map[string]int{},
		},
		{
			name:     "repeated characters",
			input:    "aaaaaa",
			expected: map[string]int{"a": 6},
		},
		{
			name:     "string with spaces and special chars",
			input:    "a b!@#",
			expected: map[string]int{"a": 1, " ": 1, "b": 1, "!": 1, "@": 1, "#": 1},
		},
		{
			name:     "unicode characters",
			input:    "café",
			expected: map[string]int{"c": 1, "a": 1, "f": 1, "é": 1},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := StringCounter(tt.input)
			if !reflect.DeepEqual(result, tt.expected) {
				t.Errorf("StringCounter(%q) = %v, want %v", tt.input, result, tt.expected)
			}
		})
	}
}

func TestUniqueElems(t *testing.T) {
	tests := []struct {
		name     string
		input    []string
		expected []string
	}{
		{
			name:     "no duplicates",
			input:    []string{"a", "b", "c"},
			expected: []string{"a", "b", "c"},
		},
		{
			name:     "with duplicates",
			input:    []string{"a", "b", "a", "c", "b"},
			expected: []string{"a", "b", "c"},
		},
		{
			name:     "empty slice",
			input:    []string{},
			expected: []string{},
		},
		{
			name:     "all same elements",
			input:    []string{"x", "x", "x"},
			expected: []string{"x"},
		},
		{
			name:     "single element",
			input:    []string{"single"},
			expected: []string{"single"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := UniqueElems(tt.input)
			if !reflect.DeepEqual(result, tt.expected) {
				t.Errorf("UniqueElems(%v) = %v, want %v", tt.input, result, tt.expected)
			}
		})
	}
}

func TestSumFloatS(t *testing.T) {
	tests := []struct {
		name     string
		input    []float64
		expected float64
	}{
		{
			name:     "positive numbers",
			input:    []float64{1.0, 2.0, 3.0},
			expected: 6.0,
		},
		{
			name:     "mixed positive and negative",
			input:    []float64{1.5, -2.5, 3.0},
			expected: 2.0,
		},
		{
			name:     "empty slice",
			input:    []float64{},
			expected: 0.0,
		},
		{
			name:     "single element",
			input:    []float64{42.5},
			expected: 42.5,
		},
		{
			name:     "zeros",
			input:    []float64{0.0, 0.0, 0.0},
			expected: 0.0,
		},
		{
			name:     "very small numbers",
			input:    []float64{1e-10, 2e-10, 3e-10},
			expected: 6e-10,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := SumFloatS(tt.input)
			if math.Abs(result-tt.expected) > 1e-15 {
				t.Errorf("SumFloatS(%v) = %v, want %v", tt.input, result, tt.expected)
			}
		})
	}
}

func TestRoundFloat(t *testing.T) {
	tests := []struct {
		name      string
		value     float64
		precision uint
		expected  float64
	}{
		{
			name:      "round to 2 decimal places",
			value:     3.14159,
			precision: 2,
			expected:  3.14,
		},
		{
			name:      "round to 0 decimal places",
			value:     3.7,
			precision: 0,
			expected:  4.0,
		},
		{
			name:      "already rounded",
			value:     5.0,
			precision: 2,
			expected:  5.0,
		},
		{
			name:      "negative number",
			value:     -2.567,
			precision: 1,
			expected:  -2.6,
		},
		{
			name:      "high precision",
			value:     1.23456789,
			precision: 5,
			expected:  1.23457,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := RoundFloat(tt.value, tt.precision)
			if math.Abs(result-tt.expected) > 1e-10 {
				t.Errorf("RoundFloat(%v, %d) = %v, want %v", tt.value, tt.precision, result, tt.expected)
			}
		})
	}
}

func TestRoundFloatS(t *testing.T) {
	tests := []struct {
		name      string
		input     []float64
		precision uint
		expected  []float64
	}{
		{
			name:      "round slice to 2 decimals",
			input:     []float64{3.14159, 2.71828, 1.41421},
			precision: 2,
			expected:  []float64{3.14, 2.72, 1.41},
		},
		{
			name:      "empty slice",
			input:     []float64{},
			precision: 2,
			expected:  []float64{},
		},
		{
			name:      "single element",
			input:     []float64{3.14159},
			precision: 1,
			expected:  []float64{3.1},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := RoundFloatS(tt.input, tt.precision)
			if len(result) != len(tt.expected) {
				t.Errorf("RoundFloatS(%v, %d) length = %d, want %d", tt.input, tt.precision, len(result), len(tt.expected))
				return
			}
			for i, val := range result {
				if math.Abs(val-tt.expected[i]) > 1e-10 {
					t.Errorf("RoundFloatS(%v, %d)[%d] = %v, want %v", tt.input, tt.precision, i, val, tt.expected[i])
				}
			}
		})
	}
}

func TestNewRational(t *testing.T) {
	tests := []struct {
		name     string
		input    float64
		checkNum int64
		checkDen int64
	}{
		{
			name:     "simple fraction",
			input:    0.5,
			checkNum: 1,
			checkDen: 2,
		},
		{
			name:     "integer",
			input:    3.0,
			checkNum: 3,
			checkDen: 1,
		},
		{
			name:     "negative fraction",
			input:    -0.25,
			checkNum: -1,
			checkDen: 4,
		},
		{
			name:     "zero",
			input:    0.0,
			checkNum: 0,
			checkDen: 1,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := NewRational(tt.input)
			result.Simplify() // Ensure it's in simplest form

			if result.Num.Int64() != tt.checkNum || result.Den.Int64() != tt.checkDen {
				t.Errorf("NewRational(%v) = %d/%d, want %d/%d",
					tt.input, result.Num.Int64(), result.Den.Int64(), tt.checkNum, tt.checkDen)
			}
		})
	}
}

func TestNewRationalSpecialCases(t *testing.T) {
	// Test infinity
	result := NewRational(math.Inf(1))
	if result.Num.Int64() != 0 || result.Den.Int64() != 1 {
		t.Errorf("NewRational(+Inf) = %d/%d, want 0/1", result.Num.Int64(), result.Den.Int64())
	}

	// Test NaN
	result = NewRational(math.NaN())
	if result.Num.Int64() != 0 || result.Den.Int64() != 1 {
		t.Errorf("NewRational(NaN) = %d/%d, want 0/1", result.Num.Int64(), result.Den.Int64())
	}
}

func TestNewRationalWithLimit(t *testing.T) {
	tests := []struct {
		name           string
		input          float64
		maxDenominator int64
		expectSimple   bool
	}{
		{
			name:           "simple fraction within limit",
			input:          0.5,
			maxDenominator: 10,
			expectSimple:   true,
		},
		{
			name:           "pi approximation",
			input:          math.Pi,
			maxDenominator: 100,
			expectSimple:   true,
		},
		{
			name:           "zero max denominator uses default",
			input:          0.333333,
			maxDenominator: 0,
			expectSimple:   true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := NewRationalWithLimit(tt.input, tt.maxDenominator)
			if result == nil {
				t.Errorf("NewRationalWithLimit(%v, %d) returned nil", tt.input, tt.maxDenominator)
			}
			if result.Den.Int64() <= 0 {
				t.Errorf("NewRationalWithLimit(%v, %d) has non-positive denominator", tt.input, tt.maxDenominator)
			}
		})
	}
}

func TestRationalSimplify(t *testing.T) {
	tests := []struct {
		name      string
		num       int64
		den       int64
		expectNum int64
		expectDen int64
	}{
		{
			name:      "already simplified",
			num:       1,
			den:       2,
			expectNum: 1,
			expectDen: 2,
		},
		{
			name:      "needs simplification",
			num:       6,
			den:       8,
			expectNum: 3,
			expectDen: 4,
		},
		{
			name:      "negative denominator",
			num:       1,
			den:       -2,
			expectNum: -1,
			expectDen: 2,
		},
		{
			name:      "both negative",
			num:       -4,
			den:       -6,
			expectNum: 2,
			expectDen: 3,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			r := &Rational{
				Num: big.NewInt(tt.num),
				Den: big.NewInt(tt.den),
			}
			r.Simplify()

			if r.Num.Int64() != tt.expectNum || r.Den.Int64() != tt.expectDen {
				t.Errorf("Simplify %d/%d = %d/%d, want %d/%d",
					tt.num, tt.den, r.Num.Int64(), r.Den.Int64(), tt.expectNum, tt.expectDen)
			}
		})
	}
}

func TestFindLCMSlice(t *testing.T) {
	tests := []struct {
		name     string
		input    []int64
		expected int64
	}{
		{
			name:     "simple case",
			input:    []int64{2, 3, 4},
			expected: 12,
		},
		{
			name:     "with common factors",
			input:    []int64{6, 8, 12},
			expected: 24,
		},
		{
			name:     "single element",
			input:    []int64{5},
			expected: 5,
		},
		{
			name:     "empty slice",
			input:    []int64{},
			expected: 1,
		},
		{
			name:     "with zero",
			input:    []int64{0, 5, 10},
			expected: 0,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			bigInts := make([]*big.Int, len(tt.input))
			for i, val := range tt.input {
				bigInts[i] = big.NewInt(val)
			}

			result := FindLCMSlice(bigInts)
			if result.Int64() != tt.expected {
				t.Errorf("FindLCMSlice(%v) = %d, want %d", tt.input, result.Int64(), tt.expected)
			}
		})
	}
}

func TestFindGCDSlice(t *testing.T) {
	tests := []struct {
		name     string
		input    []int64
		expected int64
	}{
		{
			name:     "simple case",
			input:    []int64{12, 18, 24},
			expected: 6,
		},
		{
			name:     "coprime numbers",
			input:    []int64{7, 11, 13},
			expected: 1,
		},
		{
			name:     "single element",
			input:    []int64{42},
			expected: 42,
		},
		{
			name:     "empty slice",
			input:    []int64{},
			expected: 1,
		},
		{
			name:     "with zero",
			input:    []int64{0, 15},
			expected: 15,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			bigInts := make([]*big.Int, len(tt.input))
			for i, val := range tt.input {
				bigInts[i] = big.NewInt(val)
			}

			result := FindGCDSlice(bigInts)
			if result.Int64() != tt.expected {
				t.Errorf("FindGCDSlice(%v) = %d, want %d", tt.input, result.Int64(), tt.expected)
			}
		})
	}
}

func TestSymmetricDifference(t *testing.T) {
	tests := []struct {
		name     string
		slice1   []string
		slice2   []string
		expected []string
	}{
		{
			name:     "no overlap",
			slice1:   []string{"a", "b"},
			slice2:   []string{"c", "d"},
			expected: []string{"a", "b", "c", "d"},
		},
		{
			name:     "partial overlap",
			slice1:   []string{"a", "b", "c"},
			slice2:   []string{"b", "c", "d"},
			expected: []string{"a", "d"},
		},
		{
			name:     "complete overlap",
			slice1:   []string{"a", "b"},
			slice2:   []string{"a", "b"},
			expected: []string{},
		},
		{
			name:     "empty slices",
			slice1:   []string{},
			slice2:   []string{},
			expected: []string{},
		},
		{
			name:     "one empty slice",
			slice1:   []string{"a", "b"},
			slice2:   []string{},
			expected: []string{"a", "b"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := SymmetricDifference(tt.slice1, tt.slice2)

			// Convert to maps for comparison since order might vary
			resultMap := make(map[string]bool)
			expectedMap := make(map[string]bool)

			for _, v := range result {
				resultMap[v] = true
			}
			for _, v := range tt.expected {
				expectedMap[v] = true
			}

			if !reflect.DeepEqual(resultMap, expectedMap) {
				t.Errorf("SymmetricDifference(%v, %v) = %v, want %v",
					tt.slice1, tt.slice2, result, tt.expected)
			}
		})
	}
}

func TestReplaceNthOccurrence(t *testing.T) {
	tests := []struct {
		name     string
		s        string
		old      string
		new      string
		n        int
		expected string
	}{
		{
			name:     "replace first occurrence",
			s:        "hello world hello",
			old:      "hello",
			new:      "hi",
			n:        1,
			expected: "hi world hello",
		},
		{
			name:     "replace second occurrence",
			s:        "hello world hello",
			old:      "hello",
			new:      "hi",
			n:        2,
			expected: "hello world hi",
		},
		{
			name:     "n greater than occurrences",
			s:        "hello world",
			old:      "hello",
			new:      "hi",
			n:        2,
			expected: "hello world",
		},
		{
			name:     "n is zero",
			s:        "hello world",
			old:      "hello",
			new:      "hi",
			n:        0,
			expected: "hello world",
		},
		{
			name:     "n is negative",
			s:        "hello world",
			old:      "hello",
			new:      "hi",
			n:        -1,
			expected: "hello world",
		},
		{
			name:     "empty old string",
			s:        "hello world",
			old:      "",
			new:      "hi",
			n:        1,
			expected: "hello world",
		},
		{
			name:     "old string not found",
			s:        "hello world",
			old:      "foo",
			new:      "bar",
			n:        1,
			expected: "hello world",
		},
		{
			name:     "empty string",
			s:        "",
			old:      "hello",
			new:      "hi",
			n:        1,
			expected: "",
		},
		{

			name:     "overlapping patterns",
			s:        "aaaa",
			old:      "aa",
			new:      "b",
			n:        2,
			expected: "aab",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := ReplaceNthOccurrence(tt.s, tt.old, tt.new, tt.n)
			if result != tt.expected {
				t.Errorf("ReplaceNthOccurrence(%q, %q, %q, %d) = %q, want %q",
					tt.s, tt.old, tt.new, tt.n, result, tt.expected)
			}
		})
	}
}
