package chemreaction

import (
	"context"
	"slices"
	"testing"
)

func TestCombParallelMatchesSequential(t *testing.T) {
	tests := []struct {
		reaction string
		maxCoef  uint
	}{
		{reaction: "H2+O2=H2O", maxCoef: 4},
		{reaction: "Bi2O3+TiO2=Bi4Ti3O12", maxCoef: 4},
		{reaction: "C2H6O+C2H6O=C+H2+O2", maxCoef: 4},
		{reaction: "KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O", maxCoef: 6},
		{reaction: "KMnO4+HCl=MnCl2+Cl2+H2O+KCl", maxCoef: 6},
	}
	for _, tt := range tests {
		t.Run(tt.reaction, func(t *testing.T) {
			reac, err := NewChemicalReaction(tt.reaction)
			if err != nil {
				t.Fatal(err)
			}
			bal, err := reac.Balancer()
			if err != nil {
				t.Fatal(err)
			}

			sequential := bal.bAlgos.combinatorialSequential(context.Background(), tt.maxCoef)
			parallel := bal.bAlgos.combinatorial(context.Background(), tt.maxCoef)
			if !slices.Equal(sequential, parallel) {
				t.Errorf("parallel = %v, sequential = %v", parallel, sequential)
			}
		})
	}
}
