package utils

import (
	"math"
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

func TestNewSimpleFraction(t *testing.T) {
	tests := []struct {
		name           string
		input          float64
		maxDenominator int64
		expectedNum    int64
		expectedDen    int64
	}{
		{
			name:           "simple fraction",
			input:          0.5,
			maxDenominator: 100,
			expectedNum:    1,
			expectedDen:    2,
		},
		{
			name:           "integer",
			input:          3.0,
			maxDenominator: 100,
			expectedNum:    3,
			expectedDen:    1,
		},
		{
			name:           "negative fraction",
			input:          -0.25,
			maxDenominator: 100,
			expectedNum:    -1,
			expectedDen:    4,
		},
		{
			name:           "zero",
			input:          0.0,
			maxDenominator: 100,
			expectedNum:    0,
			expectedDen:    1,
		},
		{
			name:           "one third approximation",
			input:          0.333333,
			maxDenominator: 10,
			expectedNum:    1,
			expectedDen:    3,
		},
		{
			name:           "pi approximation",
			input:          math.Pi,
			maxDenominator: 10,
			expectedNum:    22,
			expectedDen:    7,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := NewSimpleFraction(tt.input, tt.maxDenominator)
			if result.Num != tt.expectedNum || result.Den != tt.expectedDen {
				t.Errorf("NewSimpleFraction(%v, %d) = %d/%d, want %d/%d",
					tt.input, tt.maxDenominator, result.Num, result.Den, tt.expectedNum, tt.expectedDen)
			}
		})
	}
}

func TestNewSimpleFractionSpecialCases(t *testing.T) {
	result := NewSimpleFraction(math.Inf(1), 100)
	if result.Num != 0 || result.Den != 1 {
		t.Errorf("NewSimpleFraction(+Inf, 100) = %d/%d, want 0/1", result.Num, result.Den)
	}

	result = NewSimpleFraction(math.NaN(), 100)
	if result.Num != 0 || result.Den != 1 {
		t.Errorf("NewSimpleFraction(NaN, 100) = %d/%d, want 0/1", result.Num, result.Den)
	}

	result = NewSimpleFraction(math.Inf(-1), 100)
	if result.Num != 0 || result.Den != 1 {
		t.Errorf("NewSimpleFraction(-Inf, 100) = %d/%d, want 0/1", result.Num, result.Den)
	}
}

func TestGcdInt64(t *testing.T) {
	tests := []struct {
		name     string
		a        int64
		b        int64
		expected int64
	}{
		{
			name:     "simple case",
			a:        12,
			b:        8,
			expected: 4,
		},
		{
			name:     "coprime numbers",
			a:        7,
			b:        11,
			expected: 1,
		},
		{
			name:     "one is zero",
			a:        0,
			b:        5,
			expected: 5,
		},
		{
			name:     "negative numbers",
			a:        -12,
			b:        8,
			expected: 4,
		},
		{
			name:     "both negative",
			a:        -12,
			b:        -8,
			expected: 4,
		},
		{
			name:     "same numbers",
			a:        15,
			b:        15,
			expected: 15,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := gcdInt64(tt.a, tt.b)
			if result != tt.expected {
				t.Errorf("gcdInt64(%d, %d) = %d, want %d", tt.a, tt.b, result, tt.expected)
			}
		})
	}
}

func TestLcmInt64(t *testing.T) {
	tests := []struct {
		name     string
		a        int64
		b        int64
		expected int64
	}{
		{
			name:     "simple case",
			a:        4,
			b:        6,
			expected: 12,
		},
		{
			name:     "coprime numbers",
			a:        7,
			b:        11,
			expected: 77,
		},
		{
			name:     "one is zero",
			a:        0,
			b:        5,
			expected: 0,
		},
		{
			name:     "negative numbers",
			a:        -4,
			b:        6,
			expected: 12,
		},
		{
			name:     "both negative",
			a:        -4,
			b:        -6,
			expected: 12,
		},
		{
			name:     "same numbers",
			a:        15,
			b:        15,
			expected: 15,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := lcmInt64(tt.a, tt.b)
			if result != tt.expected {
				t.Errorf("lcmInt64(%d, %d) = %d, want %d", tt.a, tt.b, result, tt.expected)
			}
		})
	}
}

func TestFindLCMSliceInt64(t *testing.T) {
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
		{
			name:     "overflow case",
			input:    []int64{1e10, 1e10},
			expected: -1,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := FindLCMSliceInt64(tt.input)
			if result != tt.expected {
				t.Errorf("FindLCMSliceInt64(%v) = %d, want %d", tt.input, result, tt.expected)
			}
		})
	}
}

func TestFindGCDSliceInt64(t *testing.T) {
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
		{
			name:     "negative numbers",
			input:    []int64{-12, 18, 24},
			expected: 6,
		},
		{
			name:     "all negative",
			input:    []int64{-12, -18, -24},
			expected: 6,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := FindGCDSliceInt64(tt.input)
			if result != tt.expected {
				t.Errorf("FindGCDSliceInt64(%v) = %d, want %d", tt.input, result, tt.expected)
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
