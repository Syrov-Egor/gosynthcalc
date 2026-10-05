package chemformula

import (
	"slices"
	"testing"
)

func TestChemicalFormulaOutput(t *testing.T) {
	formulaStr := "H2SO4"
	form, _ := NewChemicalFormula(formulaStr)
	got := form.Output().String()
	expected := `formula: H2SO4
parsed formula: ['H': 2 'S': 1 'O': 4]
molar mass: 98.0720
mass percent: ['H': 2.0556 'S': 32.6903 'O': 65.2541]
atomic percent: ['H': 28.5714 'S': 14.2857 'O': 57.1429]
oxide percent: ['H2O': 18.3692 'SO3': 81.6308]`
	if got != expected {
		t.Errorf("Output() expected %s, got %s", expected, got)
	}
}

func TestChemicalFormulaOxidePercent_defaultThenCustom(t *testing.T) {
	form, err := NewChemicalFormula("BaFeO4")
	if err != nil {
		t.Fatalf("NewChemicalFormula() error: %v", err)
	}

	defaultExpected := []Atom{
		{Label: "BaO", Amount: 65.75731389},
		{Label: "Fe2O3", Amount: 34.24268611},
	}
	customExpected := []Atom{
		{Label: "BaO", Amount: 66.51800627},
		{Label: "Fe3O4", Amount: 33.48199373},
	}

	gotDefault, err := form.OxidePercent()
	if err != nil {
		t.Fatalf("OxidePercent() error: %v", err)
	}
	if !slices.Equal(gotDefault, defaultExpected) {
		t.Errorf("OxidePercent() = %v, expected %v", gotDefault, defaultExpected)
	}

	gotCustom, err := form.OxidePercent("Fe3O4")
	if err != nil {
		t.Fatalf("OxidePercent(\"Fe3O4\") error: %v", err)
	}
	if !slices.Equal(gotCustom, customExpected) {
		t.Errorf("OxidePercent(\"Fe3O4\") = %v, expected %v", gotCustom, customExpected)
	}

	gotDefaultAgain, err := form.OxidePercent()
	if err != nil {
		t.Fatalf("OxidePercent() error: %v", err)
	}
	if !slices.Equal(gotDefaultAgain, defaultExpected) {
		t.Errorf("OxidePercent() after custom call = %v, expected %v", gotDefaultAgain, defaultExpected)
	}
}

func TestChemicalFormulaOxidePercent_customThenDefault(t *testing.T) {
	form, err := NewChemicalFormula("BaFeO4")
	if err != nil {
		t.Fatalf("NewChemicalFormula() error: %v", err)
	}

	defaultExpected := []Atom{
		{Label: "BaO", Amount: 65.75731389},
		{Label: "Fe2O3", Amount: 34.24268611},
	}
	customExpected := []Atom{
		{Label: "BaO", Amount: 66.51800627},
		{Label: "Fe3O4", Amount: 33.48199373},
	}

	gotCustom, err := form.OxidePercent("Fe3O4")
	if err != nil {
		t.Fatalf("OxidePercent(\"Fe3O4\") error: %v", err)
	}
	if !slices.Equal(gotCustom, customExpected) {
		t.Errorf("OxidePercent(\"Fe3O4\") = %v, expected %v", gotCustom, customExpected)
	}

	gotDefault, err := form.OxidePercent()
	if err != nil {
		t.Fatalf("OxidePercent() error: %v", err)
	}
	if !slices.Equal(gotDefault, defaultExpected) {
		t.Errorf("OxidePercent() after custom call = %v, expected %v", gotDefault, defaultExpected)
	}

	gotOutput, err := form.OxidePercent()
	if err != nil {
		t.Fatalf("OxidePercent() error: %v", err)
	}
	if !slices.Equal(gotOutput, defaultExpected) {
		t.Errorf("OxidePercent() after custom call = %v, expected %v", gotOutput, defaultExpected)
	}
	if out := form.Output().OxidePercent; !slices.Equal(out, defaultExpected) {
		t.Errorf("Output().OxidePercent = %v, expected %v", out, defaultExpected)
	}
}
