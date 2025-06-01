package chemformula

import (
	"testing"
)

func TestChemicalFormulaOutput(t *testing.T) {
	formulaStr := "H2SO4"
	form, _ := NewChemicalFormula(formulaStr)
	got := form.Output().String()
	expected := `formula: H2SO4
parsed formula: ['H': 2 'S': 1 'O': 4]
molar mass: 98.072
mass percent: ['H': 2.0556 'S': 32.6903 'O': 65.2541]
atomic percent: ['H': 28.5714 'S': 14.2857 'O': 57.1429]
oxide percent: ['H2O': 18.3692 'SO3': 81.6308]`
	if got != expected {
		t.Errorf("Output() expected %s, got %s", expected, got)
	}
}
