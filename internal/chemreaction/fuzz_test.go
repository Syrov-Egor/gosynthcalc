package chemreaction

import (
	"strings"
	"testing"
)

func FuzzBalanceReaction(f *testing.F) {
	seeds := []string{
		"H2+O2=H2O",
		"KMnO4+HCl=MnCl2+Cl2+H2O+KCl",
		"KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O",
		"Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2",
		"H2+O2=H2O+Na",
		"FeS2+HNO3=Fe2(SO4)3+NO+H2SO4",
		"PtNH3BrNO3+CuNH3KNO3+BeCO3=C44H50O14.98+Cu(NO3)2+PtO3+Br1.99NO2+K1.97O+BeO+HNO3",
	}
	for _, seed := range seeds {
		f.Add(seed)
	}

	f.Fuzz(func(t *testing.T, reaction string) {
		if len(reaction) > 512 || strings.Count(reaction, "+") > 20 {
			t.Skip()
		}
		reac, err := NewChemicalReaction(reaction)
		if err != nil {
			return
		}
		_, _ = reac.Coefficients()
		_, _ = reac.CoefficientsExact()
		_ = reac.IsBalanced()
	})
}
