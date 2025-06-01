package chemreaction

import (
	"testing"
)

func TestReactionValidator_emptyReaction(t *testing.T) {
	tests := []struct {
		name     string
		reaction string
		expected bool
	}{
		{
			name:     "empty string",
			reaction: "",
			expected: true,
		},
		{
			name:     "non-empty string",
			reaction: "H2+O2=H2O",
			expected: false,
		},
		{
			name:     "whitespace only",
			reaction: "   ",
			expected: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := reactionValidator{reaction: tt.reaction}
			result := v.emptyReaction()
			if result != tt.expected {
				t.Errorf("emptyFormula() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestReactionValidator_invalidCharacters(t *testing.T) {
	tests := []struct {
		name     string
		reaction string
		expected []string
	}{
		{
			name:     "valid reaction with no invalid characters",
			reaction: "H2+O2=H2O",
			expected: []string{},
		},
		{
			name:     "reaction with invalid character",
			reaction: "H2+O2=H2カO",
			expected: []string{"カ"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := reactionValidator{reaction: tt.reaction}
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

func TestReactionValidator_noRPSeparator(t *testing.T) {
	tests := []struct {
		name     string
		reaction string
		expected bool
	}{
		{
			name:     "has separator",
			reaction: "H2+O2=H2O",
			expected: false,
		},
		{
			name:     "no separator",
			reaction: "H2+O2 H2O",
			expected: true,
		},
		{
			name:     "different separator",
			reaction: "H2+O2->H2O",
			expected: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := reactionValidator{reaction: tt.reaction}
			d, _ := newReactionDecomposer(tt.reaction)
			result := v.noRPSeparator(*d)
			if result != tt.expected {
				t.Errorf("noRPSeparator() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestReactionValidator_noReacSeparator(t *testing.T) {
	tests := []struct {
		name     string
		reaction string
		expected bool
	}{
		{
			name:     "has reactant separator",
			reaction: "H2 + O2 -> H2O",
			expected: false,
		},
		{
			name:     "no reactant separator",
			reaction: "H2 O2 -> H2O",
			expected: true,
		},
		{
			name:     "multiple reactant separators",
			reaction: "H2 + O2 + N2 -> H2O",
			expected: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := reactionValidator{reaction: tt.reaction}
			result := v.noReacSeparator()
			if result != tt.expected {
				t.Errorf("noReacSeparator() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestReactionValidator_validate(t *testing.T) {
	tests := []struct {
		name     string
		reaction string
		errorous bool
	}{
		{
			name:     "empty string",
			reaction: "",
			errorous: true,
		},
		{
			name:     "invalid characters",
			reaction: "K2CO3+HCl=H2CO3+KClкалий",
			errorous: true,
		},
		{
			name:     "no separator between r and p",
			reaction: "K2CO3+HCl+H2CO3+KCl",
			errorous: true,
		},
		{
			name:     "no separator between compounds",
			reaction: "K2CO 3HCl = H2CO3 KCl",
			errorous: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := reactionValidator{reaction: tt.reaction}
			_, err := v.validate()
			if err == nil {
				t.Errorf("Expected error, got nil")
			}
		})
	}
}
