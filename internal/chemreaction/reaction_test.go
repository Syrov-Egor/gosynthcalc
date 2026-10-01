package chemreaction

import (
	"slices"
	"testing"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

func TestChemicalReacutionOutput(t *testing.T) {
	reactionStr := "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
	reac, _ := NewChemicalReaction(reactionStr)
	got, _ := reac.Output()
	expected := `initial reaction: Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2
reaction matrix:
⎡2  0  3  1  1  3  0  0⎤
⎢3  0  4  1  0  0  1  2⎥
⎣0  1  0  0  0  1  1  1⎦
mode: Balance
formulas: [Fe2O3 C Fe3O4 FeO Fe Fe3C CO CO2]
coefficients: [1954 1854 518 1093 1096 55 901 898]
coefficients normalized: [3.77220077 3.57915058 1 2.11003861 2.11583012 0.10617761 1.73938224 1.73359073]
algorithm: general pseudoinverse
is balanced: true
final reaction: 1954Fe2O3+1854C=518Fe3O4+1093FeO+1096Fe+55Fe3C+901CO+898CO2
final reaction normalized: 3.77220077Fe2O3+3.57915058C=Fe3O4+2.11003861FeO+2.11583012Fe+0.10617761Fe3C+1.73938224CO+1.73359073CO2
molar masses: [159.6870 12.0110 231.5310 71.8440 55.8450 179.5460 28.0100 44.0090]
target: Fe3O4
masses: [2.6017 0.1857 1.0000 0.6547 0.5103 0.0823 0.2104 0.3295]
Fe2O3  M = 159.6870  g/mol  m = 2.6017  g
C      M = 12.0110   g/mol  m = 0.1857  g
Fe3O4  M = 231.5310  g/mol  m = 1.0000  g
FeO    M = 71.8440   g/mol  m = 0.6547  g
Fe     M = 55.8450   g/mol  m = 0.5103  g
Fe3C   M = 179.5460  g/mol  m = 0.0823  g
CO     M = 28.0100   g/mol  m = 0.2104  g
CO2    M = 44.0090   g/mol  m = 0.3295  g`
	if got.String() != expected {
		t.Errorf("Output() expected '%s', got '%s'", expected, got)
	}
}

func TestChemicalReaction_forceMode(t *testing.T) {
	reactionStr := "Cr2(SO4)3+Br2+NaOH=NaBr+Na2CrO4+Na2SO4+H2O"
	reacOpts := ReacOptions{
		Rmode:      Force,
		Target:     0,
		TargerMass: 1.0,
		Intify:     true,
		Precision:  8,
		Tolerance:  1e-8,
	}
	reac, _ := NewChemicalReaction(reactionStr, reacOpts)
	got := reac.decomposer.initCoefs
	expected := []float64{1, 1, 1, 1, 1, 1, 1}
	if !slices.Equal(got, expected) {
		t.Errorf("initCoefs expected '%v', got '%v'", expected, got)
	}
}

func TestChemicalReaction_checkModeRight(t *testing.T) {
	reactionStr := "Cr2(SO4)3+3Br2+16NaOH=6NaBr+2Na2CrO4+3Na2SO4+8H2O"
	reacOpts := ReacOptions{
		Rmode:      Check,
		Target:     0,
		TargerMass: 1.0,
		Intify:     true,
		Precision:  8,
		Tolerance:  1e-8,
	}
	reac, _ := NewChemicalReaction(reactionStr, reacOpts)
	got, _ := reac.Coefficients()
	expected := []float64{1, 3, 16, 6, 2, 3, 8}
	if !slices.Equal(got.Result, expected) {
		t.Errorf("initCoefs expected '%v', got '%v'", expected, got)
	}
}

func TestChemicalReaction_checkModeWrong(t *testing.T) {
	reactionStr := "Cr2(SO4)3+Br2+NaOH=NaBr+Na2CrO4+Na2SO4+H2O"
	reacOpts := ReacOptions{
		Rmode:      Check,
		Target:     0,
		TargerMass: 1.0,
		Intify:     true,
		Precision:  8,
		Tolerance:  1e-8,
	}
	reac, _ := NewChemicalReaction(reactionStr, reacOpts)
	_, err := reac.Coefficients()
	if err.Error() != "reaction is not balanced" {
		t.Errorf("this test should give error %s, got %s instead", "Reaction is not balanced", err)
	}
}

func TestChemicalReaction_countValidationLeft(t *testing.T) {
	reactionStr := "Rb2CO3+La2O3+Nb2O5=RbLaNb2O7"
	reac, _ := NewChemicalReaction(reactionStr)
	_, err := reac.Coefficients()
	expected := "cannot balance this reaction, because element(s) [C] are only in one part of the reaction"
	if err.Error() != expected {
		t.Errorf("this test should give error %s, got %s instead",
			expected,
			err)
	}
}

func TestChemicalReaction_countValidationRight(t *testing.T) {
	reactionStr := "Rb2CO3+La2O3+Nb2O5=RbLaNb2O7+CO2+Nd"
	reac, _ := NewChemicalReaction(reactionStr)
	_, err := reac.Coefficients()
	expected := "cannot balance this reaction, because element(s) [Nd] are only in one part of the reaction"
	if err.Error() != expected {
		t.Errorf("this test should give error %s, got %s instead",
			expected,
			err)
	}
}

func TestChemicalReaction_countValidationBoth(t *testing.T) {
	reactionStr := "Rb2CO3+La2O3+Nb2O5=RbLaNb2O7+Nd"
	reac, _ := NewChemicalReaction(reactionStr)
	_, err := reac.Coefficients()
	expected := "cannot balance this reaction, because element(s) [C Nd] are only in one part of the reaction"
	if err.Error() != expected {
		t.Errorf("this test should give error %s, got %s instead",
			expected,
			err)
	}
}

func TestChemicalReaction_setCoefficientsWronglen(t *testing.T) {
	reactionStr := "Cr2(SO4)3+Br2+NaOH=NaBr+Na2CrO4+Na2SO4+H2O"
	reac, _ := NewChemicalReaction(reactionStr)
	coefs := []float64{2, 5, 6, 1, 2, 4}
	err := reac.SetCoefficients(coefs)
	expected := "length of coefficient slice should be 7, got 6"
	if err.Error() != expected {
		t.Errorf("this test should give error %s, got %s instead",
			expected,
			err)
	}
}

func TestChemicalReaction_setCoefficientsNegative(t *testing.T) {
	reactionStr := "Cr2(SO4)3+Br2+NaOH=NaBr+Na2CrO4+Na2SO4+H2O"
	reac, _ := NewChemicalReaction(reactionStr)
	coefs := []float64{2, 5, 6, 1, 2, 4, -2}
	err := reac.SetCoefficients(coefs)
	expected := "input coefficient -2.000000 at position 6 is <= 0"
	if err.Error() != expected {
		t.Errorf("this test should give error %s, got %s instead",
			expected,
			err)
	}
}

func TestChemicalReaction_setCoefficientsRight(t *testing.T) {
	reactionStr := "Cr2(SO4)3+Br2+NaOH=NaBr+Na2CrO4+Na2SO4+H2O"
	reac, _ := NewChemicalReaction(reactionStr)
	coefs := []float64{2, 5, 6, 1, 2, 4, 2}
	reac.SetCoefficients(coefs)
	reac_coefs, _ := reac.Coefficients()
	expected := []float64{2, 5, 6, 1, 2, 4, 2}
	if !slices.Equal(reac_coefs.Result, expected) {
		t.Errorf("this test should give %v, got %v instead",
			expected,
			reac_coefs.Result)
	}
}

func TestChemicalReaction_setCoefficientsInvalidatesCaches(t *testing.T) {
	reactionStr := "H2+O2=H2O"
	reac, err := NewChemicalReaction(reactionStr)
	if err != nil {
		t.Fatalf("NewChemicalReaction() error: %v", err)
	}

	if _, err := reac.Output(); err != nil {
		t.Fatalf("Output() error: %v", err)
	}

	newCoefs := []float64{6, 3, 4}
	if err := reac.SetCoefficients(newCoefs); err != nil {
		t.Fatalf("SetCoefficients() error: %v", err)
	}

	newCoefs[0] = 999

	coefs, err := reac.Coefficients()
	if err != nil {
		t.Fatalf("Coefficients() error: %v", err)
	}
	if !slices.Equal(coefs.Result, []float64{6, 3, 4}) {
		t.Errorf("Coefficients() = %v, want [6 3 4]", coefs.Result)
	}
	if coefs.Method != "User" {
		t.Errorf("Method = %q, want \"User\"", coefs.Method)
	}

	expectedNorm := []float64{1.5, 0.75, 1}
	norm, err := reac.NormCoefficients()
	if err != nil {
		t.Fatalf("NormCoefficients() error: %v", err)
	}
	if !slices.Equal(norm, expectedNorm) {
		t.Errorf("NormCoefficients() = %v, want %v", norm, expectedNorm)
	}

	fin, err := reac.FinalReaction()
	if err != nil {
		t.Fatalf("FinalReaction() error: %v", err)
	}
	if fin != "6H2+3O2=4H2O" {
		t.Errorf("FinalReaction() = %q, want %q", fin, "6H2+3O2=4H2O")
	}

	finNorm, err := reac.FinalReactionNorm()
	if err != nil {
		t.Fatalf("FinalReactionNorm() error: %v", err)
	}
	if finNorm != "1.5H2+0.75O2=H2O" {
		t.Errorf("FinalReactionNorm() = %q, want %q", finNorm, "1.5H2+0.75O2=H2O")
	}

	molars, err := reac.MolarMasses()
	if err != nil {
		t.Fatalf("MolarMasses() error: %v", err)
	}
	target := 2
	nu := 1.0 / molars[target]
	expectedMasses := make([]float64, len(molars))
	for i, molar := range molars {
		expectedMasses[i] = utils.RoundFloat(molar*nu*expectedNorm[i], 8)
	}
	masses, err := reac.Masses()
	if err != nil {
		t.Fatalf("Masses() error: %v", err)
	}
	if !slices.Equal(masses, expectedMasses) {
		t.Errorf("Masses() = %v, want %v", masses, expectedMasses)
	}

	out, err := reac.Output()
	if err != nil {
		t.Fatalf("Output() error: %v", err)
	}
	if !slices.Equal(out.Coefficients, []float64{6, 3, 4}) {
		t.Errorf("Output().Coefficients = %v, want [6 3 4]", out.Coefficients)
	}
	if !slices.Equal(out.NormCoefficients, expectedNorm) {
		t.Errorf("Output().NormCoefficients = %v, want %v", out.NormCoefficients, expectedNorm)
	}
	if out.FinalReaction != "6H2+3O2=4H2O" {
		t.Errorf("Output().FinalReaction = %q, want %q", out.FinalReaction, "6H2+3O2=4H2O")
	}
	if out.FinalReactionNorm != "1.5H2+0.75O2=H2O" {
		t.Errorf("Output().FinalReactionNorm = %q, want %q", out.FinalReactionNorm, "1.5H2+0.75O2=H2O")
	}
	if !slices.Equal(out.Masses, expectedMasses) {
		t.Errorf("Output().Masses = %v, want %v", out.Masses, expectedMasses)
	}
}
