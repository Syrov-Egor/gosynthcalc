package chemreaction

import (
	"slices"
	"strings"
	"testing"
)

func sanitizedReactionValidator(reaction string) reactionValidator {
	return reactionValidator{reaction: sanitize(reaction)}
}

func checkValidatorResult(t *testing.T, decomp *reactionDecomposer, err error, wantErr bool, errorContains string) {
	t.Helper()
	if !wantErr {
		if err != nil {
			t.Errorf("validate() unexpected error = %v", err)
			return
		}
		if decomp == nil {
			t.Errorf("validate() returned nil decomposer and nil error")
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

func TestReactionValidator_emptyReaction(t *testing.T) {
	tests := []struct {
		name          string
		reaction      string
		wantErr       bool
		errorContains string
	}{
		{
			name:          "empty string",
			reaction:      "",
			wantErr:       true,
			errorContains: "empty reaction string",
		},
		{
			name:     "non-empty string",
			reaction: "H2+O2=H2O",
			wantErr:  false,
		},
		{
			name:          "whitespace only",
			reaction:      "   ",
			wantErr:       true,
			errorContains: "empty reaction string",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := sanitizedReactionValidator(tt.reaction)
			decomp, err := v.validate()
			checkValidatorResult(t, decomp, err, tt.wantErr, tt.errorContains)
		})
	}
}

func TestReactionValidator_invalidCharacters(t *testing.T) {
	tests := []struct {
		name          string
		reaction      string
		wantErr       bool
		errorContains string
	}{
		{
			name:     "valid reaction with no invalid characters",
			reaction: "H2+O2=H2O",
			wantErr:  false,
		},
		{
			name:          "reaction with invalid character",
			reaction:      "H2+O2=H2カO",
			wantErr:       true,
			errorContains: "invalid character(s) カ",
		},
		{
			name:          "reaction with multiple invalid characters",
			reaction:      "H2+O2=H2カO@",
			wantErr:       true,
			errorContains: "invalid character(s) カ@",
		},
		{
			name:          "raw whitespace is rejected before sanitization",
			reaction:      "H2 + O2 = H2O",
			wantErr:       true,
			errorContains: "invalid character(s)",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := reactionValidator{reaction: tt.reaction}
			decomp, err := v.validate()
			checkValidatorResult(t, decomp, err, tt.wantErr, tt.errorContains)
		})
	}
}

func TestReactionValidator_noRPSeparator(t *testing.T) {
	tests := []struct {
		name          string
		reaction      string
		wantErr       bool
		errorContains string
	}{
		{
			name:     "has separator",
			reaction: "H2+O2=H2O",
			wantErr:  false,
		},
		{
			name:          "no separator",
			reaction:      "H2+O2 H2O",
			wantErr:       true,
			errorContains: "no separator between reactants and products",
		},
		{
			name:     "different separator",
			reaction: "H2+O2->H2O",
			wantErr:  false,
		},
		{
			name:          "empty reactant side",
			reaction:      "=H2O",
			wantErr:       true,
			errorContains: "no separator between reactants and products",
		},
		{
			name:          "empty product side",
			reaction:      "H2=",
			wantErr:       true,
			errorContains: "no separator between reactants and products",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := sanitizedReactionValidator(tt.reaction)
			decomp, err := v.validate()
			checkValidatorResult(t, decomp, err, tt.wantErr, tt.errorContains)
		})
	}
}

func TestReactionValidator_noReacSeparator(t *testing.T) {
	tests := []struct {
		name          string
		reaction      string
		wantErr       bool
		errorContains string
	}{
		{
			name:     "has reactant separator",
			reaction: "H2 + O2 -> H2O",
			wantErr:  false,
		},
		{
			name:          "no reactant separator",
			reaction:      "H2 O2 -> H2O",
			wantErr:       true,
			errorContains: "no separators between compounds",
		},
		{
			name:     "multiple reactant separators",
			reaction: "H2 + O2 + N2 -> H2O",
			wantErr:  false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := sanitizedReactionValidator(tt.reaction)
			decomp, err := v.validate()
			checkValidatorResult(t, decomp, err, tt.wantErr, tt.errorContains)
		})
	}
}

func TestReactionValidator_validate(t *testing.T) {
	tests := []struct {
		name          string
		reaction      string
		wantErr       bool
		errorContains string
	}{
		{
			name:          "empty string",
			reaction:      "",
			wantErr:       true,
			errorContains: "empty reaction string",
		},
		{
			name:          "invalid characters",
			reaction:      "K2CO3+HCl=H2CO3+KClкалий",
			wantErr:       true,
			errorContains: "invalid character(s) калий",
		},
		{
			name:          "no separator between r and p",
			reaction:      "K2CO3+HCl+H2CO3+KCl",
			wantErr:       true,
			errorContains: "no separator between reactants and products",
		},
		{
			name:          "no separator between compounds",
			reaction:      "K2CO 3HCl = H2CO3 KCl",
			wantErr:       true,
			errorContains: "no separators between compounds",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := sanitizedReactionValidator(tt.reaction)
			decomp, err := v.validate()
			checkValidatorResult(t, decomp, err, tt.wantErr, tt.errorContains)
		})
	}
}

func TestSanitize(t *testing.T) {
	tests := []struct {
		name     string
		reaction string
		expected string
	}{
		{
			name:     "empty string",
			reaction: "",
			expected: "",
		},
		{
			name:     "whitespace only",
			reaction: "   ",
			expected: "",
		},
		{
			name:     "no changes needed",
			reaction: "H2+O2=H2O",
			expected: "H2+O2=H2O",
		},
		{
			name:     "spaces around separators and coefficients",
			reaction: " 2 H2 + O2 -> 2 H2O ",
			expected: "2H2+O2=2H2O",
		},
		{
			name:     "double equals separator",
			reaction: "H2==O2",
			expected: "H2=O2",
		},
		{
			name:     "reversible arrow separator",
			reaction: "H2+O2<->H2O",
			expected: "H2+O2=H2O",
		},
		{
			name:     "angle brackets separator",
			reaction: "H2+O2<>H2O",
			expected: "H2+O2=H2O",
		},
		{
			name:     "greater than separator",
			reaction: "H2+O2>H2O",
			expected: "H2+O2=H2O",
		},
		{
			name:     "unicode arrow separator",
			reaction: "H2+O2→H2O",
			expected: "H2+O2=H2O",
		},
		{
			name:     "unicode equilibrium separator",
			reaction: "H2+O2⇄H2O",
			expected: "H2+O2=H2O",
		},
		{
			name:     "adduct symbols preserved",
			reaction: "CuSO4·5H2O+NaOH=Cu(OH)2+Na2SO4",
			expected: "CuSO4·5H2O+NaOH=Cu(OH)2+Na2SO4",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := sanitize(tt.reaction)
			if result != tt.expected {
				t.Errorf("sanitize() = %q, expected %q", result, tt.expected)
			}
		})
	}
}

func TestReactionValidator_allReactionSeparators(t *testing.T) {
	for _, sep := range reactionSymbols.reactionSeparators {
		t.Run(sep, func(t *testing.T) {
			v := sanitizedReactionValidator("H2+O2" + sep + "H2O")
			decomp, err := v.validate()
			checkValidatorResult(t, decomp, err, false, "")
			if decomp == nil {
				return
			}
			if !slices.Equal(decomp.reactants, []string{"H2", "O2"}) {
				t.Errorf("reactants = %v, expected %v", decomp.reactants, []string{"H2", "O2"})
			}
			if !slices.Equal(decomp.products, []string{"H2O"}) {
				t.Errorf("products = %v, expected %v", decomp.products, []string{"H2O"})
			}
		})
	}
}

func TestReactionValidator_decomposition(t *testing.T) {
	tests := []struct {
		name          string
		reaction      string
		separatorPos  int
		initCoefs     []float64
		compounds     []string
		wantErr       bool
		errorContains string
	}{
		{
			name:         "integer coefficients",
			reaction:     "2H2+O2=2H2O",
			separatorPos: 2,
			initCoefs:    []float64{2, 1, 2},
			compounds:    []string{"H2", "O2", "H2O"},
		},
		{
			name:         "decimal coefficients",
			reaction:     "0.5N2+1.5O2=NO",
			separatorPos: 2,
			initCoefs:    []float64{0.5, 1.5, 1},
			compounds:    []string{"N2", "O2", "NO"},
		},
		{
			name:         "coefficient on product side",
			reaction:     "H2+Cl2=2HCl",
			separatorPos: 2,
			initCoefs:    []float64{1, 1, 2},
			compounds:    []string{"H2", "Cl2", "HCl"},
		},
		{
			name:          "adjacent reactant separators",
			reaction:      "H2++O2=H2O",
			wantErr:       true,
			errorContains: "two adjacent +",
		},
		{
			name:          "trailing reactant separator",
			reaction:      "H2+=H2O",
			wantErr:       true,
			errorContains: "two adjacent +",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := sanitizedReactionValidator(tt.reaction)
			decomp, err := v.validate()
			checkValidatorResult(t, decomp, err, tt.wantErr, tt.errorContains)
			if tt.wantErr || decomp == nil {
				return
			}
			if decomp.separatorPos != tt.separatorPos {
				t.Errorf("separatorPos = %v, expected %v", decomp.separatorPos, tt.separatorPos)
			}
			if !slices.Equal(decomp.initCoefs, tt.initCoefs) {
				t.Errorf("initCoefs = %v, expected %v", decomp.initCoefs, tt.initCoefs)
			}
			if !slices.Equal(decomp.compounds, tt.compounds) {
				t.Errorf("compounds = %v, expected %v", decomp.compounds, tt.compounds)
			}
		})
	}
}

func TestReactionValidator_adductSymbols(t *testing.T) {
	tests := []struct {
		name     string
		reaction string
	}{
		{
			name:     "middle dot hydrate",
			reaction: "CuSO4·5H2O+2NaOH=Cu(OH)2+Na2SO4+5H2O",
		},
		{
			name:     "bullet hydrate",
			reaction: "CuSO4•5H2O+2NaOH=Cu(OH)2+Na2SO4+5H2O",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := sanitizedReactionValidator(tt.reaction)
			decomp, err := v.validate()
			checkValidatorResult(t, decomp, err, false, "")
		})
	}
}
