package chemformula

import (
	"strings"
	"testing"
)

func TestFormulaValidator_emptyFormula(t *testing.T) {
	tests := []struct {
		name     string
		formula  string
		expected bool
	}{
		{
			name:     "empty string",
			formula:  "",
			expected: true,
		},
		{
			name:     "non-empty string",
			formula:  "H2O",
			expected: false,
		},
		{
			name:     "whitespace only",
			formula:  "   ",
			expected: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{formula: tt.formula}
			result := v.emptyFormula()
			if result != tt.expected {
				t.Errorf("emptyFormula() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestFormulaValidator_noLetters(t *testing.T) {
	tests := []struct {
		name     string
		formula  string
		expected bool
	}{
		{
			name:     "empty formula",
			formula:  "",
			expected: true,
		},
		{
			name:     "just brackets",
			formula:  "[]",
			expected: true,
		},
		{
			name:     "non-empty string",
			formula:  "H2O",
			expected: false,
		},
		{
			name:     "just numbers",
			formula:  "222",
			expected: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{formula: tt.formula}
			result := v.noLetters()
			if result != tt.expected {
				t.Errorf("noLetters() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestFormulaValidator_invalidCharacters(t *testing.T) {
	tests := []struct {
		name     string
		formula  string
		expected []string
	}{
		{
			name:     "valid formula with no invalid characters",
			formula:  "H2O",
			expected: []string{},
		},
		{
			name:     "formula with special characters",
			formula:  "H2O@#$",
			expected: []string{"@", "#", "$"},
		},
		{
			name:     "formula with valid brackets and adducts",
			formula:  "Ca(OH)2·H2O",
			expected: []string{},
		},
		{
			name:     "formula with invalid punctuation",
			formula:  "H2O,NH3",
			expected: []string{","},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{formula: tt.formula}
			result := v.invalidCharacters()
			if len(result) != len(tt.expected) {
				t.Errorf("invalidCharacters() = %v, expected %v", result, tt.expected)
				return
			}
			for i, char := range result {
				if i >= len(tt.expected) || char != tt.expected[i] {
					t.Errorf("invalidCharacters() = %v, expected %v", result, tt.expected)
					break
				}
			}
		})
	}
}

func TestFormulaValidator_invalidAtoms(t *testing.T) {
	tests := []struct {
		name     string
		formula  string
		expected []string
	}{
		{
			name:     "valid atoms only",
			formula:  "H2O",
			expected: []string{},
		},
		{
			name:     "invalid single atom",
			formula:  "Xy2O",
			expected: []string{"Xy"},
		},
		{
			name:     "multiple invalid atoms",
			formula:  "XyZw3",
			expected: []string{"Xy", "Zw"},
		},
		{
			name:     "valid complex formula",
			formula:  "Ca(OH)2",
			expected: []string{},
		},
		{
			name:     "mixed valid and invalid",
			formula:  "NaClXy",
			expected: []string{"Xy"},
		},
		{
			name:     "test leftover 1",
			formula:  "Li(ac)*2H2O",
			expected: []string{"a", "c"},
		},
		{
			name:     "test leftover 2",
			formula:  "aLi*2H2O",
			expected: []string{"a"},
		},
		{
			name:     "test invalid atoms and leftovers",
			formula:  "ALk*2H2O",
			expected: []string{"A", "Lk"},
		},
		{
			name:     "test overlapping atoms 1",
			formula:  "OsPoPO3",
			expected: []string{},
		},
		{
			name:     "test overlapping atoms 2",
			formula:  "[Ru(C10H8N2)3]Cl2*6H2O",
			expected: []string{},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{formula: tt.formula}
			result := v.invalidAtoms()

			resultMap := make(map[string]bool)
			for _, atom := range result {
				resultMap[atom] = true
			}

			expectedMap := make(map[string]bool)
			for _, atom := range tt.expected {
				expectedMap[atom] = true
			}

			if len(resultMap) != len(expectedMap) {
				t.Errorf("invalidAtoms() = %v, expected %v", result, tt.expected)
				return
			}

			for atom := range expectedMap {
				if !resultMap[atom] {
					t.Errorf("invalidAtoms() = %v, expected %v", result, tt.expected)
					break
				}
			}
		})
	}
}

func TestFormulaValidator_bracketsBalance(t *testing.T) {
	tests := []struct {
		name     string
		formula  string
		expected bool
	}{
		{
			name:     "balanced parentheses",
			formula:  "Ca(OH)2",
			expected: true,
		},
		{
			name:     "unbalanced parentheses - missing closing",
			formula:  "Ca(OH2",
			expected: false,
		},
		{
			name:     "unbalanced parentheses - missing opening",
			formula:  "CaOH)2",
			expected: false,
		},
		{
			name:     "balanced square brackets",
			formula:  "K3[Fe(CN)6]",
			expected: true,
		},
		{
			name:     "unbalanced square brackets",
			formula:  "K3[Fe(CN)6",
			expected: false,
		},
		{
			name:     "balanced curly brackets",
			formula:  "Cu{NH3}4",
			expected: true,
		},
		{
			name:     "mixed balanced brackets",
			formula:  "K3[Fe(CN)6]{H2O}",
			expected: true,
		},
		{
			name:     "mixed unbalanced brackets",
			formula:  "K3[Fe(CN)6]{H2O",
			expected: false,
		},
		{
			name:     "no brackets",
			formula:  "NaCl",
			expected: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{formula: tt.formula}
			result := v.bracketsBalance()
			if result != tt.expected {
				t.Errorf("bracketsBalance() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestFormulaValidator_numOfAdducts(t *testing.T) {
	tests := []struct {
		name     string
		formula  string
		expected int
	}{
		{
			name:     "no adducts",
			formula:  "H2O",
			expected: 0,
		},
		{
			name:     "single dot adduct",
			formula:  "CaCl2·6H2O",
			expected: 1,
		},
		{
			name:     "single bullet adduct",
			formula:  "CaCl2•6H2O",
			expected: 1,
		},
		{
			name:     "multiple adducts",
			formula:  "CaCl2·6H2O·NaCl",
			expected: 2,
		},
		{
			name:     "mixed adduct symbols",
			formula:  "CaCl2·6H2O•NaCl",
			expected: 2,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{formula: tt.formula}
			result := v.numOfAdducts()
			if result != tt.expected {
				t.Errorf("numOfAdducts() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestFormulaValidator_validate(t *testing.T) {
	tests := []struct {
		name          string
		formula       string
		expectError   bool
		errorContains string
	}{
		{
			name:        "valid simple formula",
			formula:     "H2O",
			expectError: false,
		},
		{
			name:        "valid complex formula",
			formula:     "Ca(OH)2",
			expectError: false,
		},
		{
			name:          "empty formula",
			formula:       "",
			expectError:   true,
			errorContains: "Empty formula string",
		},
		{
			name:          "invalid characters",
			formula:       "H2O@",
			expectError:   true,
			errorContains: "invalid character(s)",
		},
		{
			name:          "invalid atoms",
			formula:       "XyO2",
			expectError:   true,
			errorContains: "invalid atom(s)",
		},
		{
			name:          "unbalanced brackets",
			formula:       "Ca(OH2",
			expectError:   true,
			errorContains: "not balanced",
		},
		{
			name:          "multiple adducts",
			formula:       "CaCl2·6H2O·NaCl",
			expectError:   true,
			errorContains: "more than 1 adduct",
		},
		{
			name:        "single adduct valid",
			formula:     "CaCl2·6H2O",
			expectError: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{formula: tt.formula}
			err := v.validate()

			if tt.expectError {
				if err == nil {
					t.Errorf("validate() expected error, got nil")
					return
				}
				if tt.errorContains != "" && !strings.Contains(err.Error(), tt.errorContains) {
					t.Errorf("validate() error = %v, expected to contain %v", err.Error(), tt.errorContains)
				}
			} else {
				if err != nil {
					t.Errorf("validate() unexpected error = %v", err)
				}
			}
		})
	}
}

func BenchmarkFormulaValidator_validate(b *testing.B) {
	testCases := []string{
		"H2O",
		"Ca(OH)2",
		"K3[Fe(CN)6]",
		"CaCl2·6H2O",
		"C6H12O6",
	}

	for _, formula := range testCases {
		b.Run(formula, func(b *testing.B) {
			v := formulaValidator{formula: formula}
			for b.Loop() {
				_ = v.validate()
			}
		})
	}
}
