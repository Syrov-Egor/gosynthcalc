package chemreaction

import (
	"context"
	"testing"

	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
)

func BenchmarkCombSequential(b *testing.B) {
	benchmarkComb(b, false)
}

func BenchmarkCombParallel(b *testing.B) {
	benchmarkComb(b, true)
}

func benchmarkComb(b *testing.B, parallel bool) {
	reac, err := NewChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl")
	if err != nil {
		b.Fatal(err)
	}
	bal, err := reac.Balancer()
	if err != nil {
		b.Fatal(err)
	}

	ctx := context.Background()
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var result []float64
		if parallel {
			result = bal.bAlgos.combinatorial(ctx, 16)
		} else {
			result = bal.bAlgos.combinatorialSequential(ctx, 16)
		}
		if result == nil {
			b.Fatal("Comb found no solution")
		}
	}
}

func BenchmarkAutoCorpus(b *testing.B) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		b.Fatal(err)
	}

	type reactor struct {
		reaction *ChemicalReaction
		balancer *balancer
	}
	reactors := make([]reactor, 0, len(reactions))
	for _, rd := range reactions {
		reac, err := NewChemicalReaction(rd.reaction)
		if err != nil {
			b.Fatalf("NewChemicalReaction(%q): %v", rd.reaction, err)
		}
		bal, err := reac.Balancer()
		if err != nil {
			b.Fatalf("Balancer(%q): %v", rd.reaction, err)
		}
		// Warm the exact analysis cache so this benchmark measures routing
		// and solving rather than the one-time float-to-rational conversion.
		if _, err := bal.analysis(); err != nil {
			b.Fatalf("analysis(%q): %v", rd.reaction, err)
		}
		reactors = append(reactors, reactor{reaction: reac, balancer: bal})
	}

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for _, r := range reactors {
			if _, err := r.balancer.Auto(); err != nil {
				b.Fatalf("Auto(%q): %v", r.reaction.reaction, err)
			}
		}
	}
}

func BenchmarkExactAnalysisCorpus(b *testing.B) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		b.Fatal(err)
	}

	type item struct {
		parsed       [][]chemformula.ExactAtom
		separatorPos int
	}
	items := make([]item, 0, len(reactions))
	for _, rd := range reactions {
		reac, err := NewChemicalReaction(rd.reaction)
		if err != nil {
			b.Fatalf("NewChemicalReaction(%q): %v", rd.reaction, err)
		}
		parsed, err := reac.ExactParsedFormulas()
		if err != nil {
			b.Fatalf("ExactParsedFormulas(%q): %v", rd.reaction, err)
		}
		items = append(items, item{parsed: parsed, separatorPos: reac.decomposer.separatorPos})
	}

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for _, it := range items {
			exact, err := createExactReacMatrix(it.parsed, it.separatorPos)
			if err != nil {
				b.Fatalf("createExactReacMatrix: %v", err)
			}
			if _, err := analyze(exact); err != nil {
				b.Fatalf("analyze: %v", err)
			}
		}
	}
}

func BenchmarkExactSimplexCorpus(b *testing.B) {
	reactions, err := parseReactionsCSV("testing_reactions.csv")
	if err != nil {
		b.Fatal(err)
	}

	balancers := make([]*balancer, 0, len(reactions))
	for _, rd := range reactions {
		reac, err := NewChemicalReaction(rd.reaction)
		if err != nil {
			b.Fatalf("NewChemicalReaction(%q): %v", rd.reaction, err)
		}
		bal, err := reac.Balancer()
		if err != nil {
			b.Fatalf("Balancer(%q): %v", rd.reaction, err)
		}
		if _, err := bal.analysis(); err != nil {
			b.Fatalf("analysis(%q): %v", rd.reaction, err)
		}
		balancers = append(balancers, bal)
	}

	ctx := context.Background()
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for _, bal := range balancers {
			if _, err := bal.UniversalExact(ctx); err != nil {
				b.Fatalf("UniversalExact: %v", err)
			}
		}
	}
}
