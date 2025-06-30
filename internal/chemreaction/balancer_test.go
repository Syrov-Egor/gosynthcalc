package chemreaction

import (
	"context"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"log"
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

func parseReactionsCSV(filename string) ([]reactionData, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()

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
