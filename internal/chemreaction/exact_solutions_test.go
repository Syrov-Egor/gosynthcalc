package chemreaction

import (
	"encoding/csv"
	"fmt"
	"math/big"
	"os"
	"strings"
	"testing"
)

type exactSolutionRow struct {
	reaction string
	method   string
	coefs    []*big.Int
}

func parseExactSolutionsCSV(filename string) ([]exactSolutionRow, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer func() {
		_ = file.Close()
	}()

	records, err := csv.NewReader(file).ReadAll()
	if err != nil {
		return nil, err
	}
	if len(records) < 2 {
		return nil, fmt.Errorf("exact-solutions CSV is empty")
	}

	rows := make([]exactSolutionRow, 0, len(records)-1)
	for i, record := range records[1:] {
		if len(record) != 3 {
			return nil, fmt.Errorf("row %d has %d fields, want 3", i+2, len(record))
		}
		coefs, err := parseBigIntList(record[2])
		if err != nil {
			return nil, fmt.Errorf("row %d: %w", i+2, err)
		}
		rows = append(rows, exactSolutionRow{
			reaction: strings.TrimSpace(record[0]),
			method:   strings.TrimSpace(record[1]),
			coefs:    coefs,
		})
	}
	return rows, nil
}

func parseBigIntList(s string) ([]*big.Int, error) {
	s = strings.TrimSpace(s)
	if !strings.HasPrefix(s, "[") || !strings.HasSuffix(s, "]") {
		return nil, fmt.Errorf("invalid coefficient list %q", s)
	}
	s = strings.TrimSpace(s[1 : len(s)-1])
	if s == "" {
		return []*big.Int{}, nil
	}
	parts := strings.Split(s, ",")
	out := make([]*big.Int, len(parts))
	for i, part := range parts {
		part = strings.TrimSpace(part)
		v, ok := new(big.Int).SetString(part, 10)
		if !ok {
			return nil, fmt.Errorf("invalid integer %q", part)
		}
		out[i] = v
	}
	return out, nil
}

func TestExactSolutionsGoldenCorpus(t *testing.T) {
	rows, err := parseExactSolutionsCSV("testing_exact_solutions.csv")
	if err != nil {
		t.Fatal(err)
	}
	if len(rows) != 113 {
		t.Fatalf("golden corpus has %d rows, want 113", len(rows))
	}

	for i, row := range rows {
		reac, err := NewChemicalReaction(row.reaction)
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		got, err := reac.CoefficientsExact()
		if err != nil {
			t.Fatalf("row %d: CoefficientsExact: %v", i+1, err)
		}
		if got.Method != row.method {
			t.Errorf("row %d: method = %q, want %q", i+1, got.Method, row.method)
		}
		if len(got.Result) != len(row.coefs) {
			t.Errorf("row %d: coefficient count = %d, want %d",
				i+1, len(got.Result), len(row.coefs))
			continue
		}
		for j := range row.coefs {
			if got.Result[j].Cmp(row.coefs[j]) != 0 {
				t.Errorf("row %d coefficient %d = %s, want %s",
					i+1, j, got.Result[j], row.coefs[j])
			}
		}

		bal, err := reac.Balancer()
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		A, _, err := bal.exactMatrixAndAnalysis()
		if err != nil {
			t.Fatalf("row %d: %v", i+1, err)
		}
		requirePrimitivePositive(t, got.Result, A)
	}
}
