package chemformula

import (
	"math/big"
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

func TestNewChemicalFormula_parsesAndCaches(t *testing.T) {
	tests := []struct {
		formula  string
		expected []Atom
	}{
		{"H2O", []Atom{{Label: "H", Amount: 2}, {Label: "O", Amount: 1}}},
		{"K3[Fe(CN)6]", []Atom{
			{Label: "K", Amount: 3}, {Label: "Fe", Amount: 1},
			{Label: "C", Amount: 6}, {Label: "N", Amount: 6},
		}},
		{"(K0,6Na0,4)2[S]O4", []Atom{
			{Label: "K", Amount: 1.2}, {Label: "Na", Amount: 0.8},
			{Label: "S", Amount: 1}, {Label: "O", Amount: 4},
		}},
		{" {Cu[O H]2} · 0,5 H2 O ", []Atom{
			{Label: "Cu", Amount: 1}, {Label: "O", Amount: 2.5}, {Label: "H", Amount: 3},
		}},
		{"H.5O", []Atom{{Label: "H", Amount: 0.5}, {Label: "O", Amount: 1}}},
		{"H0O", []Atom{{Label: "H", Amount: 0}, {Label: "O", Amount: 1}}},
	}
	for _, tt := range tests {
		t.Run(tt.formula, func(t *testing.T) {
			formula, err := NewChemicalFormula(tt.formula, 6)
			if err != nil {
				t.Fatalf("NewChemicalFormula() unexpected error: %v", err)
			}
			if !slices.Equal(formula.parsedFormula, tt.expected) {
				t.Errorf("constructor did not cache parsed result: got %v, expected %v", formula.parsedFormula, tt.expected)
			}
			first, second := formula.ParsedFormula(), formula.ParsedFormula()
			if !slices.Equal(first, tt.expected) {
				t.Errorf("ParsedFormula() = %v, expected %v", first, tt.expected)
			}
			if len(first) == 0 || len(second) == 0 || &first[0] != &second[0] {
				t.Error("ParsedFormula() did not reuse the cached result")
			}
			if formula.Formula() != tt.formula || formula.sanFormula != sanitize(tt.formula) || formula.precision != 6 {
				t.Error("constructor did not preserve formula, normalization, or precision")
			}
		})
	}
}

func TestChemicalFormula_ExactParsedFormula(t *testing.T) {
	tests := []struct {
		formula string
		want    []ExactAtom
	}{
		{
			formula: "(Ga0.9Sc0.1)3",
			want: []ExactAtom{
				{Label: "Ga", Amount: big.NewRat(27, 10)},
				{Label: "Sc", Amount: big.NewRat(3, 10)},
			},
		},
		{
			formula: "(H.5O0.25)1.5",
			want: []ExactAtom{
				{Label: "H", Amount: big.NewRat(3, 4)},
				{Label: "O", Amount: big.NewRat(3, 8)},
			},
		},
		{
			formula: "CuSO4*.5(H2O)2NaCl",
			want: []ExactAtom{
				{Label: "Cu", Amount: big.NewRat(1, 1)},
				{Label: "S", Amount: big.NewRat(1, 1)},
				{Label: "O", Amount: big.NewRat(5, 1)},
				{Label: "H", Amount: big.NewRat(2, 1)},
				{Label: "Na", Amount: big.NewRat(1, 2)},
				{Label: "Cl", Amount: big.NewRat(1, 2)},
			},
		},
		{
			formula: "H0O",
			want: []ExactAtom{
				{Label: "H", Amount: new(big.Rat)},
				{Label: "O", Amount: big.NewRat(1, 1)},
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.formula, func(t *testing.T) {
			formula, err := NewChemicalFormula(tt.formula)
			if err != nil {
				t.Fatalf("NewChemicalFormula() error: %v", err)
			}
			got := formula.ExactParsedFormula()
			if len(got) != len(tt.want) {
				t.Fatalf("ExactParsedFormula() length = %d, want %d", len(got), len(tt.want))
			}
			for i := range tt.want {
				if got[i].Label != tt.want[i].Label {
					t.Errorf("atom %d label = %q, want %q", i, got[i].Label, tt.want[i].Label)
				}
				if got[i].Amount == nil || got[i].Amount.Cmp(tt.want[i].Amount) != 0 {
					t.Errorf("atom %d amount = %v, want %v", i, got[i].Amount, tt.want[i].Amount)
				}
			}
			if !slices.EqualFunc(got, formula.ExactParsedFormula(), func(a, b ExactAtom) bool {
				return a.Label == b.Label && a.Amount.Cmp(b.Amount) == 0
			}) {
				t.Error("ExactParsedFormula() did not reuse the cached result")
			}
		})
	}
}

func TestNewChemicalFormula_invalid(t *testing.T) {
	for _, input := range []string{
		"", "   ", "H1.2.3O", "H.O", "H.", "H1.", "H1,2,3O", "(H)1.2.3",
		"(H).", "H*1.2.3O", "2H2O", ")H(", "H()", "H*", "*H", "H*2", "(H*O)",
		"[H)", "{H]", "([H)]", "]H[", "H[O", "H٢O", "H2@O",
	} {
		t.Run(input, func(t *testing.T) {
			formula, err := NewChemicalFormula(input)
			if err == nil {
				t.Fatal("NewChemicalFormula() expected error, got nil")
			}
			if formula != nil {
				t.Errorf("NewChemicalFormula() returned a formula on failure: %v", formula)
			}
		})
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
