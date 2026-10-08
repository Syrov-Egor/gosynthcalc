package chemreaction

import (
	"context"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"log"
	"math"
	"os"
	"slices"
	"strconv"
	"strings"
	"testing"
)

type reactionData struct {
	reaction string
	coefs    []float64
}

func TestBalancer_Inv(t *testing.T) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		log.Fatal(err)
	}
	for _, reaction := range reactions[:100] {
		t.Logf("%v", reaction.reaction)
		reac, _ := NewChemicalReaction(reaction.reaction)
		bal, _ := reac.Balancer()
		inv, _ := bal.Inv()
		if !slices.Equal(inv, reaction.coefs) {
			t.Errorf("Inv() method fault for reaction %v: expected %v, got %v",
				reaction.reaction,
				reaction.coefs,
				inv)
		}
	}
}

func TestBalancer_GPInv(t *testing.T) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		log.Fatal(err)
	}
	for _, reaction := range reactions[101:107] {
		t.Logf("%v", reaction.reaction)
		reac, _ := NewChemicalReaction(reaction.reaction)
		bal, _ := reac.Balancer()
		inv, _ := bal.GPinv()
		if !slices.Equal(inv, reaction.coefs) {
			t.Errorf("GPInv() method fault for reaction %v: expected %v, got %v",
				reaction.reaction,
				reaction.coefs,
				inv)
		}
	}
}

func TestBalancer_PPInv(t *testing.T) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		log.Fatal(err)
	}
	for _, reaction := range reactions[108:110] {
		t.Logf("%v", reaction.reaction)
		reac, _ := NewChemicalReaction(reaction.reaction)
		bal, _ := reac.Balancer()
		inv, _ := bal.PPinv()
		if !slices.Equal(inv, reaction.coefs) {
			t.Errorf("PPInv() method fault for reaction %v: expected %v, got %v",
				reaction.reaction,
				reaction.coefs,
				inv)
		}
	}
}

func TestBalancer_Comb(t *testing.T) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		log.Fatal(err)
	}
	for _, reaction := range reactions[111:113] {
		t.Logf("%v", reaction.reaction)
		reac, _ := NewChemicalReaction(reaction.reaction)
		bal, _ := reac.Balancer()
		inv, _ := bal.Comb(context.Background(), 10)
		if !slices.Equal(inv, reaction.coefs) {
			t.Errorf("Comb() method fault for reaction %v: expected %v, got %v",
				reaction.reaction,
				reaction.coefs,
				inv)
		}
	}
}

func TestBalancer_IntifyDoesNotClampAboveLimit(t *testing.T) {
	// The exact integer solution needs a coefficient above the supported
	// limit (1_000_000). It must never be clamped down to the limit, because
	// that silently unbalances the reaction and corrupts downstream masses.
	reac, err := NewChemicalReaction("H+O2=H1000001O2")
	if err != nil {
		t.Fatalf("NewChemicalReaction() error: %v", err)
	}

	coefs, err := reac.Coefficients()
	if err != nil {
		t.Fatalf("Coefficients() error: %v", err)
	}

	want := []float64{1000001, 1, 1}
	if !slices.Equal(coefs.Result, want) {
		t.Errorf("Coefficients() = %v, want %v", coefs.Result, want)
	}
	if !reac.IsBalanced() {
		t.Errorf("reaction H+O2=H1000001O2 with coefficients %v is not balanced", coefs.Result)
	}
}

func TestBalancer_IntifyAtLimit(t *testing.T) {
	// The exact integer solution at the supported limit must still be
	// integerified instead of falling back to floating point coefficients.
	reac, err := NewChemicalReaction("H+O2=H1000000O2")
	if err != nil {
		t.Fatalf("NewChemicalReaction() error: %v", err)
	}

	coefs, err := reac.Coefficients()
	if err != nil {
		t.Fatalf("Coefficients() error: %v", err)
	}

	want := []float64{1000000, 1, 1}
	if !slices.Equal(coefs.Result, want) {
		t.Errorf("Coefficients() = %v, want %v", coefs.Result, want)
	}
	if !reac.IsBalanced() {
		t.Errorf("reaction H+O2=H1000000O2 with coefficients %v is not balanced", coefs.Result)
	}
}

func TestBalancer_intifyCoefsPreservesUnrepresentableValues(t *testing.T) {
	reac, err := NewChemicalReaction("H+O2=H2O")
	if err != nil {
		t.Fatalf("NewChemicalReaction() error: %v", err)
	}
	bal, err := reac.Balancer()
	if err != nil {
		t.Fatalf("Balancer() error: %v", err)
	}

	tests := []struct {
		name  string
		coefs []float64
		want  []float64
	}{
		{
			name:  "below the limit",
			coefs: []float64{2, 1, 2},
			want:  []float64{2, 1, 2},
		},
		{
			name:  "at the limit",
			coefs: []float64{1000000, 1, 1},
			want:  []float64{1000000, 1, 1},
		},
		{
			name:  "above the limit keeps validated floats",
			coefs: []float64{1000001, 1, 1},
			want:  []float64{1000001, 1, 1},
		},
		{
			name:  "non finite value keeps validated floats",
			coefs: []float64{math.NaN(), 1, 1},
			want:  []float64{math.NaN(), 1, 1},
		},
		{
			name:  "negative coefficient keeps validated floats",
			coefs: []float64{-2.5, 1, 2},
			want:  []float64{-2.5, 1, 2},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			input := slices.Clone(tt.coefs)
			got := bal.intifyCoefs(input, bal.maxDenom)
			if !equalCoefs(got, tt.want) {
				t.Errorf("intifyCoefs(%v) = %v, want %v", tt.coefs, got, tt.want)
			}
			if !equalCoefs(input, tt.coefs) {
				t.Errorf("intifyCoefs mutated its input: got %v, want %v", input, tt.coefs)
			}
		})
	}
}

func equalCoefs(got, want []float64) bool {
	if len(got) != len(want) {
		return false
	}
	for i := range got {
		if math.IsNaN(want[i]) {
			if !math.IsNaN(got[i]) {
				return false
			}
			continue
		}
		if got[i] != want[i] {
			return false
		}
	}
	return true
}

func parseReactionsCSV(filename string) ([]reactionData, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %v", err)
	}

	reader := csv.NewReader(file)
	reader.TrimLeadingSpace = true

	records, err := reader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("error reading CSV: %v", err)
	}

	var reactions []reactionData

	for i, record := range records {
		if i == 0 {
			continue // Skip header
		}

		if len(record) == 0 || record[0] == "" {
			continue
		}

		reaction := strings.TrimSpace(record[0])
		coefsStr := strings.TrimSpace(record[1])

		coefsStr = strings.Trim(coefsStr, "\"")

		coefs, err := parseCoefficients(coefsStr)
		if err != nil {
			log.Printf("Warning: Failed to parse coefficients for reaction '%s': %v", reaction, err)
			continue
		}

		reactions = append(reactions, reactionData{
			reaction: reaction,
			coefs:    coefs,
		})
	}

	fileCloseErr := file.Close()
	if fileCloseErr != nil {
		return nil, err
	}
	return reactions, nil
}

func parseCoefficients(coefsStr string) ([]float64, error) {
	var intCoefs []int
	err := json.Unmarshal([]byte(coefsStr), &intCoefs)
	if err == nil {
		floatCoefs := make([]float64, len(intCoefs))
		for i, v := range intCoefs {
			floatCoefs[i] = float64(v)
		}
		return floatCoefs, nil
	}
	cleaned := strings.Trim(coefsStr, "[]")
	parts := strings.Split(cleaned, ",")

	var coefs []float64
	for _, part := range parts {
		part = strings.TrimSpace(part)
		if part == "" {
			continue
		}

		val, err := strconv.ParseFloat(part, 64)
		if err != nil {
			return nil, fmt.Errorf("invalid coefficient '%s': %v", part, err)
		}
		coefs = append(coefs, val)
	}

	return coefs, nil
}
