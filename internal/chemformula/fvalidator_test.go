package chemformula

import (
	"slices"
	"strings"
	"testing"
)

func sanitizedFormulaValidator(formula string) formulaValidator {
	return formulaValidator{formula: sanitize(formula)}
}

func checkValidatorErr(t *testing.T, err error, wantErr bool, errorContains string) {
	t.Helper()
	if !wantErr {
		if err != nil {
			t.Errorf("validate() unexpected error = %v", err)
		}
		return
	}
	if err == nil {
		t.Errorf("validate() expected error, got nil")
		return
	}
	if errorContains != "" && !strings.Contains(err.Error(), errorContains) {
		t.Errorf("validate() error = %v, expected to contain %v", err.Error(), errorContains)
	}
}

func TestFormulaValidator_emptyFormula(t *testing.T) {
	tests := []struct {
		name          string
		formula       string
		wantErr       bool
		errorContains string
	}{
		{
			name:          "empty string",
			formula:       "",
			wantErr:       true,
			errorContains: "Empty formula string",
		},
		{
			name:    "non-empty string",
			formula: "H2O",
			wantErr: false,
		},
		{
			name:          "whitespace only",
			formula:       "   ",
			wantErr:       true,
			errorContains: "Empty formula string",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := sanitizedFormulaValidator(tt.formula)
			checkValidatorErr(t, v.validate(), tt.wantErr, tt.errorContains)
		})
	}
}

func TestFormulaValidator_noLetters(t *testing.T) {
	tests := []struct {
		name          string
		formula       string
		wantErr       bool
		errorContains string
	}{
		{
			name:          "empty formula",
			formula:       "",
			wantErr:       true,
			errorContains: "Empty formula string",
		},
		{
			name:          "just brackets",
			formula:       "[]",
			wantErr:       true,
			errorContains: "No letters A-Z or a-z",
		},
		{
			name:    "non-empty string",
			formula: "H2O",
			wantErr: false,
		},
		{
			name:          "just numbers",
			formula:       "222",
			wantErr:       true,
			errorContains: "No letters A-Z or a-z",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := sanitizedFormulaValidator(tt.formula)
			checkValidatorErr(t, v.validate(), tt.wantErr, tt.errorContains)
		})
	}
}

func TestFormulaValidator_invalidCharacters(t *testing.T) {
	tests := []struct {
		name     string
		formula  string
		sanitize bool
		expected []rune
	}{
		{
			name:     "valid formula with no invalid characters",
			formula:  "H2O",
			expected: nil,
		},
		{
			name:     "formula with special characters",
			formula:  "H2O@#$",
			expected: []rune("@#$"),
		},
		{
			name:     "formula with valid brackets and adducts",
			formula:  "Ca(OH)2·H2O",
			sanitize: true,
			expected: nil,
		},
		{
			name:     "formula with invalid punctuation",
			formula:  "H2O,NH3",
			expected: []rune(","),
		},
		{
			name:     "raw brackets and adducts must be sanitized",
			formula:  "K3[Fe(CN)6]·H2O",
			expected: []rune("[]·"),
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			f := tt.formula
			if tt.sanitize {
				f = sanitize(f)
			}
			v := formulaValidator{formula: f}
			var result []rune
			for _, r := range f {
				if !v.isAllowed(r) {
					result = append(result, r)
				}
			}
			if !slices.Equal(result, tt.expected) {
				t.Errorf("invalid characters = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestFormulaValidator_isAllowed(t *testing.T) {
	tests := []struct {
		name     string
		r        rune
		expected bool
	}{
		{name: "uppercase letter", r: 'H', expected: true},
		{name: "lowercase letter", r: 'e', expected: true},
		{name: "digit", r: '2', expected: true},
		{name: "open parenthesis", r: '(', expected: true},
		{name: "close parenthesis", r: ')', expected: true},
		{name: "adduct symbol", r: '*', expected: true},
		{name: "decimal separator", r: '.', expected: true},
		{name: "at sign", r: '@', expected: false},
		{name: "hash sign", r: '#', expected: false},
		{name: "dollar sign", r: '$', expected: false},
		{name: "comma", r: ',', expected: false},
		{name: "space", r: ' ', expected: false},
		{name: "square bracket", r: '[', expected: false},
		{name: "close square bracket", r: ']', expected: false},
		{name: "curly bracket", r: '{', expected: false},
		{name: "close curly bracket", r: '}', expected: false},
		{name: "middle dot", r: '·', expected: false},
		{name: "bullet", r: '•', expected: false},
		{name: "non-ASCII rune", r: 'Ω', expected: false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{}
			result := v.isAllowed(tt.r)
			if result != tt.expected {
				t.Errorf("isAllowed(%q) = %v, expected %v", tt.r, result, tt.expected)
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
		{
			name:     "deuterium-like symbol is not an element",
			formula:  "D2O",
			expected: []string{"D"},
		},
		{
			name:     "lowercase leading atom",
			formula:  "hCl",
			expected: []string{"h"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := invalidAtoms(tt.formula)

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
		name          string
		formula       string
		wantErr       bool
		errorContains string
	}{
		{
			name:    "balanced parentheses",
			formula: "Ca(OH)2",
		},
		{
			name:          "unbalanced parentheses - missing closing",
			formula:       "Ca(OH2",
			wantErr:       true,
			errorContains: "not balanced",
		},
		{
			name:          "unbalanced parentheses - missing opening",
			formula:       "CaOH)2",
			wantErr:       true,
			errorContains: "not balanced",
		},
		{
			name:    "balanced square brackets",
			formula: "K3[Fe(CN)6]",
		},
		{
			name:          "unbalanced square brackets",
			formula:       "K3[Fe(CN)6",
			wantErr:       true,
			errorContains: "not balanced",
		},
		{
			name:    "balanced curly brackets",
			formula: "Cu{NH3}4",
		},
		{
			name:    "mixed balanced brackets",
			formula: "K3[Fe(CN)6]{H2O}",
		},
		{
			name:          "mixed unbalanced brackets",
			formula:       "K3[Fe(CN)6]{H2O",
			wantErr:       true,
			errorContains: "not balanced",
		},
		{
			name:    "no brackets",
			formula: "NaCl",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := sanitizedFormulaValidator(tt.formula)
			checkValidatorErr(t, v.validate(), tt.wantErr, tt.errorContains)
		})
	}
}

func TestFormulaValidator_numOfAdducts(t *testing.T) {
	tests := []struct {
		name          string
		formula       string
		wantErr       bool
		errorContains string
	}{
		{
			name:    "no adducts",
			formula: "H2O",
		},
		{
			name:    "single dot adduct",
			formula: "CaCl2·6H2O",
		},
		{
			name:    "single bullet adduct",
			formula: "CaCl2•6H2O",
		},
		{
			name:          "multiple adducts",
			formula:       "CaCl2·6H2O·NaCl",
			wantErr:       true,
			errorContains: "more than 1 adduct",
		},
		{
			name:          "mixed adduct symbols",
			formula:       "CaCl2·6H2O•NaCl",
			wantErr:       true,
			errorContains: "more than 1 adduct",
		},
		{
			name:    "raw adduct symbol",
			formula: "CaCl2*6H2O",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := sanitizedFormulaValidator(tt.formula)
			checkValidatorErr(t, v.validate(), tt.wantErr, tt.errorContains)
		})
	}
}

func TestFormulaValidator_validate(t *testing.T) {
	tests := []struct {
		name          string
		formula       string
		wantErr       bool
		errorContains string
	}{
		{
			name:    "valid simple formula",
			formula: "H2O",
		},
		{
			name:    "valid complex formula",
			formula: "Ca(OH)2",
		},
		{
			name:          "empty formula",
			formula:       "",
			wantErr:       true,
			errorContains: "Empty formula string",
		},
		{
			name:          "invalid characters",
			formula:       "H2O@",
			wantErr:       true,
			errorContains: "invalid character(s)",
		},
		{
			name:          "invalid atoms",
			formula:       "XyO2",
			wantErr:       true,
			errorContains: "invalid atom(s)",
		},
		{
			name:          "no letters",
			formula:       "222",
			wantErr:       true,
			errorContains: "No letters A-Z or a-z",
		},
		{
			name:          "unbalanced brackets",
			formula:       "Ca(OH2",
			wantErr:       true,
			errorContains: "not balanced",
		},
		{
			name:          "multiple adducts",
			formula:       "CaCl2·6H2O·NaCl",
			wantErr:       true,
			errorContains: "more than 1 adduct",
		},
		{
			name:    "single adduct valid",
			formula: "CaCl2·6H2O",
		},
		{
			name:    "valid coordination compound",
			formula: "K3[Fe(CN)6]",
		},
		{
			name:    "valid hydrate with square brackets",
			formula: "[Cu(NH3)4]SO4·5H2O",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := sanitizedFormulaValidator(tt.formula)
			checkValidatorErr(t, v.validate(), tt.wantErr, tt.errorContains)
		})
	}
}

func TestSanitize(t *testing.T) {
	tests := []struct {
		name     string
		formula  string
		expected string
	}{
		{
			name:     "empty string",
			formula:  "",
			expected: "",
		},
		{
			name:     "valid formula unchanged",
			formula:  "H2O",
			expected: "H2O",
		},
		{
			name:     "whitespace removed",
			formula:  "H 2 O",
			expected: "H2O",
		},
		{
			name:     "whitespace only removed",
			formula:  "   ",
			expected: "",
		},
		{
			name:     "square brackets to parentheses",
			formula:  "K3[Fe(CN)6]",
			expected: "K3(Fe(CN)6)",
		},
		{
			name:     "curly brackets to parentheses",
			formula:  "Cu{NH3}4",
			expected: "Cu(NH3)4",
		},
		{
			name:     "middle dot to adduct symbol",
			formula:  "Ca(OH)2·H2O",
			expected: "Ca(OH)2*H2O",
		},
		{
			name:     "bullet to adduct symbol",
			formula:  "CaCl2•6H2O",
			expected: "CaCl2*6H2O",
		},
		{
			name:     "comma to dot",
			formula:  "H2O,NH3",
			expected: "H2O.NH3",
		},
		{
			name:     "raw adduct symbol unchanged",
			formula:  "Ca(OH)2*H2O",
			expected: "Ca(OH)2*H2O",
		},
		{
			name:     "mixed normalization",
			formula:  "[Cu(NH3)4]SO4·5H2O",
			expected: "(Cu(NH3)4)SO4*5H2O",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := sanitize(tt.formula)
			if result != tt.expected {
				t.Errorf("sanitize() = %q, expected %q", result, tt.expected)
			}
		})
	}
}

func TestIsLetter(t *testing.T) {
	tests := []struct {
		name     string
		r        rune
		expected bool
	}{
		{name: "uppercase letter", r: 'H', expected: true},
		{name: "lowercase letter", r: 'e', expected: true},
		{name: "digit", r: '2', expected: false},
		{name: "space", r: ' ', expected: false},
		{name: "bracket", r: '(', expected: false},
		{name: "middle dot", r: '·', expected: false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := isLetter(tt.r)
			if result != tt.expected {
				t.Errorf("isLetter(%q) = %v, expected %v", tt.r, result, tt.expected)
			}
		})
	}
}

func TestIsValidElement(t *testing.T) {
	tests := []struct {
		name     string
		token    string
		expected bool
	}{
		{name: "single letter element", token: "H", expected: true},
		{name: "two letter element", token: "He", expected: true},
		{name: "two letter element 2", token: "Au", expected: true},
		{name: "superheavy element", token: "Og", expected: true},
		{name: "made up symbol", token: "Xy", expected: false},
		{name: "made up symbol 2", token: "Zw", expected: false},
		{name: "made up symbol 3", token: "Lk", expected: false},
		{name: "made up symbol 4", token: "Qq", expected: false},
		{name: "letter not used as symbol", token: "A", expected: false},
		{name: "letter not used as symbol 2", token: "Jl", expected: false},
		{name: "lowercase only", token: "h", expected: false},
		{name: "wrong case", token: "hE", expected: false},
		{name: "all uppercase", token: "HE", expected: false},
		{name: "symbol with digit", token: "H2", expected: false},
		{name: "three letters", token: "Uuo", expected: false},
		{name: "empty token", token: "", expected: false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := isValidElement(tt.token)
			if result != tt.expected {
				t.Errorf("isValidElement(%q) = %v, expected %v", tt.token, result, tt.expected)
			}
		})
	}
}

func TestIsValidElement_periodicTableComplete(t *testing.T) {
	for _, symbol := range periodicTableElements {
		t.Run(symbol, func(t *testing.T) {
			if !isValidElement(symbol) {
				t.Errorf("isValidElement(%q) = false, expected true for element from periodic table", symbol)
			}
		})
	}
}
