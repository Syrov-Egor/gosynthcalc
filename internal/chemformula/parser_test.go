package chemformula

import (
	"math"
	"slices"
	"strings"
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
			name:    "nested brackets",
			formula: "(K2)2Mg2((SO4)3Ho)2",
			expected: []Atom{
				{Label: "K", Amount: 4},
				{Label: "Mg", Amount: 2},
				{Label: "S", Amount: 6},
				{Label: "O", Amount: 24},
				{Label: "Ho", Amount: 2}},
		},
		{
			name:    "square bracket multiplier",
			formula: "[SO4]2",
			expected: []Atom{
				{Label: "S", Amount: 2},
				{Label: "O", Amount: 8}},
		},
		{
			name:    "fractional group and element amounts",
			formula: "(H.5O0.25)1.5",
			expected: []Atom{
				{Label: "H", Amount: 0.75},
				{Label: "O", Amount: 0.375}},
		},
		{
			name:    "fractional adduct multiplier applies to entire sequence",
			formula: "CuSO4*.5(H2O)2NaCl",
			expected: []Atom{
				{Label: "Cu", Amount: 1},
				{Label: "S", Amount: 1},
				{Label: "O", Amount: 5},
				{Label: "H", Amount: 2},
				{Label: "Na", Amount: 0.5},
				{Label: "Cl", Amount: 0.5}},
		},
		{
			name:    "repeated elements retain first appearance order",
			formula: "OH2(CH3OH)2*2H2O",
			expected: []Atom{
				{Label: "O", Amount: 5},
				{Label: "H", Amount: 14},
				{Label: "C", Amount: 2}},
		},
		{
			name:    "zero amounts are retained",
			formula: "H0(O2)0*0NaCl",
			expected: []Atom{
				{Label: "H", Amount: 0},
				{Label: "O", Amount: 0},
				{Label: "Na", Amount: 0},
				{Label: "Cl", Amount: 0}},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := NewParser(sanitize(tt.formula))
			result, err := v.parse()
			if err != nil {
				t.Fatalf("parse() unexpected error = %v", err)
			}
			if !slices.Equal(result, tt.expected) {
				t.Errorf("parse() = %v, expected %v", result, tt.expected)
			}
			if v.current.Type != TokenEOF {
				t.Errorf("parse() did not consume entire input: %v", v.current)
			}
		})
	}
}

func TestChemicalFormulaParser_invalid(t *testing.T) {
	tests := []struct {
		formula       string
		errorContains string
	}{
		{"", "Empty formula string"},
		{"H1.2.3O", "invalid number"},
		{"H.O", "invalid number"},
		{"H.", "invalid number"},
		{"H1.", "invalid number"},
		{"H..2O", "invalid number"},
		{"H2O1.2.3", "invalid number"},
		{"(H)1.2.3O", "invalid number"},
		{"(H).", "invalid number"},
		{"(H)2.", "invalid number"},
		{"H*1.2.3O", "invalid number"},
		{"H*.", "invalid number"},
		{"H*2.O", "invalid number"},
		{"H*2..O", "invalid number"},
		{"2H2O", "Unexpected number"},
		{".5H2O", "Unexpected number"},
		{"(2H)", "Unexpected number"},
		{"222", "No letters"},
		{")H(", "not balanced"},
		{"H)O(", "not balanced"},
		{"(H", "not balanced"},
		{"H)", "not balanced"},
		{"((H)", "not balanced"},
		{"(H))", "not balanced"},
		{"H()", "Empty parentheses group"},
		{"(())", "Empty parentheses group"},
		{"*H", "nonempty formula sequence"},
		{"H*", "nonempty formula sequence"},
		{"H*2", "nonempty formula sequence"},
		{"H**O", "nonempty formula sequence"},
		{"H*2*O", "nonempty formula sequence"},
		{"H*H*O", "more than 1 adduct"},
		{"(H*O)", "top level"},
		{"H*(O*H)", "top level"},
		{"H2@O", "invalid character(s)"},
		{"H-2", "invalid character(s)"},
		{"H+2", "invalid character(s)"},
		{"H\tO", "invalid character(s)"},
		{"H\nO", "invalid character(s)"},
		{"H٢", "invalid character(s)"},
		{"ΩH", "invalid character(s)"},
		{"Hé", "invalid character(s)"},
		{"hCl", "invalid atom(s)"},
		{"H2o", "invalid atom(s)"},
		{"H1e2", "invalid atom(s)"},
		{"Xy2O", "invalid atom(s)"},
		{"D2O", "invalid atom(s)"},
		{"H[O]2", "invalid character(s)"},
	}
	for _, tt := range tests {
		t.Run(tt.formula, func(t *testing.T) {
			result, err := NewParser(tt.formula).parse()
			if err == nil {
				t.Fatalf("parse() expected error, got %v", result)
			}
			if result != nil {
				t.Errorf("parse() returned partial result on failure: %v", result)
			}
			if !strings.Contains(err.Error(), tt.errorContains) {
				t.Errorf("parse() error = %v, expected to contain %q", err, tt.errorContains)
			}
		})
	}
}

func TestChemicalFormulaParser_numericOverflow(t *testing.T) {
	tooLarge := strings.Repeat("9", 309)
	large := "1" + strings.Repeat("0", 308)
	tests := []struct {
		name          string
		formula       string
		errorContains string
	}{
		{"element number", "H" + tooLarge, "invalid number"},
		{"group multiplier", "H(O)" + tooLarge, "invalid number"},
		{"adduct multiplier", "H*" + tooLarge + "O", "invalid number"},
		{"group multiplication", "(H2)" + large, "non-finite atom amount"},
		{"nested multiplication", "((H)" + large + ")2", "non-finite atom amount"},
		{"element accumulation", "H" + large + "H" + large, "non-finite atom amount"},
		{"group accumulation", "H" + large + "(H)" + large, "non-finite atom amount"},
		{"adduct accumulation", "H" + large + "*" + large + "H", "non-finite atom amount"},
		{"adduct multiplication", "O*" + large + "H2", "non-finite atom amount"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result, err := NewParser(tt.formula).parse()
			if err == nil || !strings.Contains(err.Error(), tt.errorContains) {
				t.Errorf("parse() error = %v, expected to contain %q", err, tt.errorContains)
			}
			if result != nil {
				t.Errorf("parse() returned partial result on failure: %v", result)
			}
		})
	}

	result, err := NewParser("H" + large).parse()
	if err != nil {
		t.Fatalf("parse() rejected finite amount: %v", err)
	}
	if len(result) != 1 || result[0].Amount != 1e308 {
		t.Errorf("parse() = %v, expected H with amount 1e308", result)
	}
}

func TestChemicalFormulaParser_errorPosition(t *testing.T) {
	tests := []struct {
		formula string
		token   string
		at      string
	}{
		{"H1.2.3O", "1.2.3", "position 2"},
		{"(H)1.2.3", "1.2.3", "position 4"},
		{"H*1.2.3O", "1.2.3", "position 3"},
		{"H2@O", "@", "position 3"},
		{"H٢O", "٢", "position 2"},
	}
	for _, tt := range tests {
		t.Run(tt.formula, func(t *testing.T) {
			_, err := NewParser(tt.formula).parse()
			if err == nil {
				t.Fatal("parse() expected error, got nil")
			}
			for _, part := range []string{tt.formula, tt.token, tt.at} {
				if !strings.Contains(err.Error(), part) {
					t.Errorf("parse() error = %v, expected to contain %q", err, part)
				}
			}
		})
	}
}

func TestLexerNextToken_invalidCharactersAreNotSkipped(t *testing.T) {
	lexer := NewLexer([]rune("H2@O٢*."))
	expected := []Token{
		{Type: TokenElement, Value: "H", Position: 0},
		{Type: TokenNumber, Value: "2", Position: 1},
		{Type: TokenInvalid, Value: "@", Position: 2},
		{Type: TokenElement, Value: "O", Position: 3},
		{Type: TokenInvalid, Value: "٢", Position: 4},
		{Type: TokenAdduct, Value: "*", Position: 5},
		{Type: TokenNumber, Value: ".", Position: 6},
		{Type: TokenEOF, Position: 7},
	}
	for _, want := range expected {
		if got := lexer.NextToken(); got != want {
			t.Errorf("NextToken() = %+v, expected %+v", got, want)
		}
	}
}

func FuzzChemicalFormulaParser_parse(f *testing.F) {
	for _, formula := range []string{
		"H2O", "(K0.6Na0.4)2(S)O4", "((H2O)2)3", "CuSO4*.5(H2O)2",
		"H0O", "H1.2.3O", "H.O", "H.", ")H(", "H*", "(H*O)", "H٢", "H2@O",
	} {
		f.Add(formula)
	}
	f.Fuzz(func(t *testing.T, formula string) {
		parser := NewParser(formula)
		atoms, err := parser.parse()
		if err != nil {
			if atoms != nil {
				t.Fatalf("parse() returned partial result on failure: %v", atoms)
			}
			return
		}
		if len(atoms) == 0 || parser.current.Type != TokenEOF {
			t.Fatalf("parse() accepted an empty or incompletely parsed formula %q", formula)
		}
		seen := make(map[string]bool)
		for _, atom := range atoms {
			if !isValidElement(atom.Label) || seen[atom.Label] {
				t.Errorf("parse() returned an invalid or duplicate atom: %v", atom)
			}
			if atom.Amount < 0 || math.IsNaN(atom.Amount) || math.IsInf(atom.Amount, 0) {
				t.Errorf("parse() returned an invalid amount: %v", atom)
			}
			seen[atom.Label] = true
		}
	})
}
