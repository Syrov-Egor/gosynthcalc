package chemreaction

import (
	"context"
	"errors"
	"io"
	"os"
	"slices"
	"testing"
)

func TestGPinvGZeroCases(t *testing.T) {
	tests := []struct {
		reaction string
		want     []float64
	}{
		{reaction: "H2+O2=H2O", want: []float64{2, 1, 2}},
		{reaction: "Bi2O3+TiO2=Bi4Ti3O12", want: []float64{2, 3, 1}},
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
			got, err := bal.GPinv()
			if err != nil {
				t.Fatalf("GPinv() error: %v", err)
			}
			if !slices.Equal(got, tt.want) {
				t.Errorf("GPinv() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestGPinvRankZero(t *testing.T) {
	m := matDenseFromFloats(t, [][]float64{
		{0, 0, 0},
		{0, 0, 0},
	})
	bal := newBalancer(m, 2, true, 8, StrategyExact, 1e-8)
	got, err := bal.GPinv()
	if err != nil {
		t.Fatalf("GPinv() error: %v", err)
	}
	want := []float64{1, 1, 1}
	if !slices.Equal(got, want) {
		t.Errorf("GPinv() = %v, want %v", got, want)
	}
}

func TestGPinvNegativeProjection(t *testing.T) {
	reac, err := NewChemicalReaction("Fe2O2H+FeO2H2=Fe3+Fe+Fe2O3H2")
	if err != nil {
		t.Fatal(err)
	}
	bal, err := reac.Balancer()
	if err != nil {
		t.Fatal(err)
	}
	_, err = bal.GPinv()
	if !errors.Is(err, ErrMethodNotPositive) {
		t.Errorf("GPinv() error = %v, want ErrMethodNotPositive", err)
	}
}

func TestPPinvGZeroCases(t *testing.T) {
	tests := []struct {
		reaction string
		want     []float64
	}{
		{reaction: "H2+O2=H2O", want: []float64{2, 1, 2}},
		{reaction: "Bi2O3+TiO2=Bi4Ti3O12", want: []float64{2, 3, 1}},
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
			got, err := bal.PPinv()
			if err != nil {
				t.Fatalf("PPinv() error: %v", err)
			}
			if !slices.Equal(got, tt.want) {
				t.Errorf("PPinv() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestPPinvRankZero(t *testing.T) {
	m := matDenseFromFloats(t, [][]float64{
		{0, 0, 0},
		{0, 0, 0},
	})
	bal := newBalancer(m, 2, true, 8, StrategyExact, 1e-8)
	got, err := bal.PPinv()
	if err != nil {
		t.Fatalf("PPinv() error: %v", err)
	}
	want := []float64{1, 1, 1}
	if !slices.Equal(got, want) {
		t.Errorf("PPinv() = %v, want %v", got, want)
	}
}

func TestPPinvNonPositive(t *testing.T) {
	reactions := []string{
		"Fe2O2H+FeO2H2=Fe3+Fe+Fe2O3H2",
		"H2+O2=H2O+Na",
	}
	for _, reaction := range reactions {
		t.Run(reaction, func(t *testing.T) {
			reac, err := NewChemicalReaction(reaction)
			if err != nil {
				t.Fatal(err)
			}
			bal, err := reac.Balancer()
			if err != nil {
				t.Fatal(err)
			}
			_, err = bal.PPinv()
			if !errors.Is(err, ErrMethodNotPositive) {
				t.Errorf("PPinv() error = %v, want ErrMethodNotPositive", err)
			}
		})
	}
}

func TestCombIsSilent(t *testing.T) {
	reac, err := NewChemicalReaction("H2+O2=H2O")
	if err != nil {
		t.Fatal(err)
	}
	bal, err := reac.Balancer()
	if err != nil {
		t.Fatal(err)
	}

	oldStdout := os.Stdout
	reader, writer, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	os.Stdout = writer
	_, combErr := bal.Comb(context.Background(), 2)
	_ = writer.Close()
	os.Stdout = oldStdout

	output, readErr := io.ReadAll(reader)
	_ = reader.Close()
	if readErr != nil {
		t.Fatal(readErr)
	}
	if len(output) != 0 {
		t.Errorf("Comb wrote to stdout: %q", output)
	}
	if combErr != nil {
		t.Errorf("Comb() error: %v", combErr)
	}
}

func TestCombCancellation(t *testing.T) {
	reac, err := NewChemicalReaction("H2+O2=H2O")
	if err != nil {
		t.Fatal(err)
	}
	bal, err := reac.Balancer()
	if err != nil {
		t.Fatal(err)
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	_, err = bal.Comb(ctx, 2)
	if !errors.Is(err, context.Canceled) {
		t.Errorf("Comb() error = %v, want context.Canceled", err)
	}
}

func TestMatrixRankRelativeTolerance(t *testing.T) {
	// A tiny but non-zero matrix must keep rank 1 under a relative
	// tolerance; an absolute 1e-8 threshold would incorrectly call it rank 0.
	m := matDenseFromFloats(t, [][]float64{{1e-12}})
	rank, _, err := matrixRank(m, 1e-8)
	if err != nil {
		t.Fatal(err)
	}
	if rank != 1 {
		t.Errorf("matrixRank(tiny matrix) = %d, want 1", rank)
	}

	zero := matDenseFromFloats(t, [][]float64{{0, 0}, {0, 0}})
	rank, _, err = matrixRank(zero, 1e-8)
	if err != nil {
		t.Fatal(err)
	}
	if rank != 0 {
		t.Errorf("matrixRank(zero matrix) = %d, want 0", rank)
	}

	if _, _, err := matrixRank(nil, 1e-8); !errors.Is(err, ErrInvalidMatrix) {
		t.Errorf("matrixRank(nil) error = %v, want ErrInvalidMatrix", err)
	}
}

func TestComputePseudoinverseRankZero(t *testing.T) {
	m := matDenseFromFloats(t, [][]float64{{1, 2}, {3, 4}})
	inverse, err := computePseudoinverse(m, 0)
	if err != nil {
		t.Fatal(err)
	}
	rows, cols := inverse.Dims()
	if rows != 2 || cols != 2 {
		t.Fatalf("zero pseudoinverse dims = %dx%d, want 2x2", rows, cols)
	}
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			if inverse.At(i, j) != 0 {
				t.Errorf("zero pseudoinverse[%d][%d] = %v, want 0", i, j, inverse.At(i, j))
			}
		}
	}

	if _, err := computePseudoinverse(m, 3); err == nil {
		t.Error("computePseudoinverse with rank above min dimension should fail")
	}
	if _, err := computePseudoinverse(nil, 0); !errors.Is(err, ErrInvalidMatrix) {
		t.Errorf("computePseudoinverse(nil) error = %v, want ErrInvalidMatrix", err)
	}
}
