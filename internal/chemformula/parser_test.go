package chemformula

import (
	"slices"
	"testing"
)

func TestChemicalFormulaParser_parse(t *testing.T) {
	tests := []struct {
		name     string
		formula  string
		expected []Atom
	}{
		{
			name:    "simple formula",
			formula: "H2O",
			expected: []Atom{
				{Label: "H", Amount: 2},
				{Label: "O", Amount: 1}},
		},
		{
			name:    "brackets with float amounts",
			formula: "(K0.6Na0.4)2[S]O4",
			expected: []Atom{
				{Label: "K", Amount: 1.2},
				{Label: "Na", Amount: 0.8},
				{Label: "S", Amount: 1},
				{Label: "O", Amount: 4}},
		},
		{
			name:    "adduct",
			formula: "(NH4)2SO4*H2O",
			expected: []Atom{
				{Label: "N", Amount: 2},
				{Label: "H", Amount: 10},
				{Label: "S", Amount: 1},
				{Label: "O", Amount: 5}},
		},
		{
			name:    "all brackets",
			formula: "{K2}2Mg2[(SO4)3Ho]2",
			expected: []Atom{
				{Label: "K", Amount: 4},
				{Label: "Mg", Amount: 2},
				{Label: "S", Amount: 6},
				{Label: "O", Amount: 24},
				{Label: "Ho", Amount: 2}},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := chemicalFormulaParser{}
			result := v.parse(tt.formula)
			if !slices.Equal(result, tt.expected) {
				t.Errorf("parse() = %v, expected %v", result, tt.expected)
			}
		})
	}
}
