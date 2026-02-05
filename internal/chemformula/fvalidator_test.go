package chemformula

import (
	"errors"
	"slices"
	"testing"
)

func Test_formulaValidator_validate_EmptyFormulaError(t *testing.T) {
	tests := []struct {
		name    string
		formula string
		wantErr bool
	}{
		{
			name:    "empty string",
			formula: "",
			wantErr: true,
		},
		{
			name:    "non-empty string",
			formula: "H2O",
			wantErr: false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{
				initialFormula: tt.formula,
				sanFormula:     tt.formula,
			}
			gotErr := v.validate()
			if tt.wantErr && gotErr == nil {
				t.Errorf("Expected error, got nil")
			}
			if !tt.wantErr && gotErr != nil {
				t.Errorf("Expected no error, got %v", gotErr)
			}
			if tt.wantErr {
				var emptyErr EmptyFormulaError
				if !errors.As(gotErr, &emptyErr) {
					t.Errorf("Expected EmptyFormulaError, got %T", gotErr)
				}
			}
		})
	}
}

func Test_formulaValidator_validate_NoLettersError(t *testing.T) {
	tests := []struct {
		name    string
		formula string
		wantErr bool
	}{
		{
			name:    "just brackets",
			formula: "[]",
			wantErr: true,
		},
		{
			name:    "just numbers",
			formula: "222",
			wantErr: true,
		},
		{
			name:    "non-empty string",
			formula: "H2O",
			wantErr: false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{
				initialFormula: tt.formula,
				sanFormula:     tt.formula,
			}
			gotErr := v.validate()
			if tt.wantErr && gotErr == nil {
				t.Errorf("Expected error, got nil")
			}
			if !tt.wantErr && gotErr != nil {
				t.Errorf("Expected no error, got %v", gotErr)
			}
			if tt.wantErr {
				var emptyErr NoLettersError
				if !errors.As(gotErr, &emptyErr) {
					t.Errorf("Expected NoLettersError, got %T", gotErr)
				}
			}
		})
	}
}

func Test_formulaValidator_validate_InvalidSymbolsError(t *testing.T) {
	tests := []struct {
		name          string
		formula       string
		wantErr       bool
		wantedSymbols []string
	}{
		{
			name:          "formula with special characters",
			formula:       "H2O@#$",
			wantErr:       true,
			wantedSymbols: []string{"@", "#", "$"},
		},
		{
			name:          "formula with cyrillic characters",
			formula:       "Hг2O",
			wantErr:       true,
			wantedSymbols: []string{"г"},
		},
		{
			name:          "formula with invalid punctuation",
			formula:       "H2O,NH3",
			wantErr:       true,
			wantedSymbols: []string{","},
		},
		{
			name:          "formula with valid brackets and adducts",
			formula:       "Ca(OH)2*H2O",
			wantErr:       false,
			wantedSymbols: []string{},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{
				initialFormula: tt.formula,
				sanFormula:     tt.formula,
			}
			gotErr := v.validate()
			if tt.wantErr && gotErr == nil {
				t.Errorf("Expected error, got nil")
			}
			if !tt.wantErr && gotErr != nil {
				t.Errorf("Expected no error, got %v", gotErr)
			}
			if tt.wantErr {
				var emptyErr InvalidSymbolsError
				if !errors.As(gotErr, &emptyErr) {
					t.Errorf("Expected InvalidSymbolsError, got %T", gotErr)
				}
				if !slices.Equal(emptyErr.symbols, tt.wantedSymbols) {
					t.Errorf("Expected %v symbols, got %v", tt.wantedSymbols, emptyErr.symbols)
				}
			}
		})
	}
}

func Test_formulaValidator_validate_InvalidAtomsError(t *testing.T) {
	tests := []struct {
		name     string
		formula  string
		wantErr  bool
		expected []string
	}{
		{
			name:     "valid atoms only",
			formula:  "H2O",
			wantErr:  false,
			expected: []string{},
		},
		{
			name:     "invalid single atom",
			formula:  "Xy2O",
			wantErr:  true,
			expected: []string{"Xy"},
		},
		{
			name:     "multiple invalid atoms",
			formula:  "XyZw3",
			wantErr:  true,
			expected: []string{"Xy", "Zw"},
		},
		{
			name:     "valid complex formula",
			formula:  "Ca(OH)2",
			wantErr:  false,
			expected: []string{},
		},
		{
			name:     "mixed valid and invalid",
			formula:  "NaClXy",
			wantErr:  true,
			expected: []string{"Xy"},
		},
		{
			name:     "test leftover 1",
			formula:  "Li(ac)*2H2O",
			wantErr:  true,
			expected: []string{"a", "c"},
		},
		{
			name:     "Abc formula",
			formula:  "Abc",
			wantErr:  true,
			expected: []string{"Ab", "c"},
		},
		{
			name:     "test leftover 2",
			formula:  "aLi*2H2O",
			wantErr:  true,
			expected: []string{"a"},
		},
		{
			name:     "test invalid atoms and leftovers",
			formula:  "ALk*2H2O",
			wantErr:  true,
			expected: []string{"A", "Lk"},
		},
		{
			name:     "test overlapping atoms 1",
			formula:  "OsPoPO3",
			wantErr:  false,
			expected: []string{},
		},
		{
			name:     "test overlapping atoms 2",
			formula:  "(Ru(C10H8N2)3)Cl2*6H2O",
			wantErr:  false,
			expected: []string{},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{
				initialFormula: tt.formula,
				sanFormula:     tt.formula,
			}
			gotErr := v.validate()
			if tt.wantErr && gotErr == nil {
				t.Errorf("Expected error, got nil")
			}
			if !tt.wantErr && gotErr != nil {
				t.Errorf("Expected no error, got %v", gotErr)
			}
			if tt.wantErr {
				var emptyErr InvalidAtomsError
				if !errors.As(gotErr, &emptyErr) {
					t.Errorf("Expected InvalidAtomsError, got %T", gotErr)
				}
				if !slices.Equal(emptyErr.atoms, tt.expected) {
					t.Errorf("Expected %v symbols, got %v", tt.expected, emptyErr.atoms)
				}
			}
		})
	}
}

func Test_formulaValidator_validate_BracketsNotBalancedError(t *testing.T) {
	tests := []struct {
		name    string
		formula string
		wantErr bool
	}{
		{
			name:    "balanced parentheses",
			formula: "Ca(OH)2",
			wantErr: false,
		},
		{
			name:    "unbalanced parentheses - missing closing",
			formula: "Ca(OH2",
			wantErr: true,
		},
		{
			name:    "unbalanced parentheses - missing opening",
			formula: "CaOH)2",
			wantErr: true,
		},
		{
			name:    "no brackets",
			formula: "NaCl",
			wantErr: false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{
				initialFormula: tt.formula,
				sanFormula:     tt.formula,
			}
			gotErr := v.validate()
			if tt.wantErr && gotErr == nil {
				t.Errorf("Expected error, got nil")
			}
			if !tt.wantErr && gotErr != nil {
				t.Errorf("Expected no error, got %v", gotErr)
			}
			if tt.wantErr {
				var emptyErr BracketsNotBalancedError
				if !errors.As(gotErr, &emptyErr) {
					t.Errorf("Expected BracketsNotBalancedError, got %T", gotErr)
				}
			}
		})
	}
}

func Test_formulaValidator_validate_MoreThanOneAdductError(t *testing.T) {
	tests := []struct {
		name    string
		formula string
		wantErr bool
	}{
		{
			name:    "no adducts",
			formula: "H2O",
			wantErr: false,
		},
		{
			name:    "single dot adduct",
			formula: "CaCl2*6H2O",
			wantErr: false,
		},
		{
			name:    "multiple adducts",
			formula: "CaCl2*6H2O*NaCl",
			wantErr: true,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := formulaValidator{
				initialFormula: tt.formula,
				sanFormula:     tt.formula,
			}
			gotErr := v.validate()
			if tt.wantErr && gotErr == nil {
				t.Errorf("Expected error, got nil")
			}
			if !tt.wantErr && gotErr != nil {
				t.Errorf("Expected no error, got %v", gotErr)
			}
			if tt.wantErr {
				var emptyErr MoreThanOneAdductError
				if !errors.As(gotErr, &emptyErr) {
					t.Errorf("Expected MoreThanOneAdductError, got %T", gotErr)
				}
			}
		})
	}
}
