package chemreaction

// Method names reported in [MethodResult.Method] and [crOutput.Algorithm].
const (
	// MethodUser marks coefficients taken from the reaction string or set
	// with [ChemicalReaction.SetCoefficients].
	MethodUser = "User"
	// MethodInverse is the exact nullity-one solver behind [balancer.Inv].
	MethodInverse = "inverse"
	// MethodExactLP is the exact rational Phase-I simplex.
	MethodExactLP = "exact LP"
	// MethodGPInv is the legacy general-pseudoinverse heuristic.
	MethodGPInv = "general pseudoinverse"
	// MethodPPInv is the legacy partial-pseudoinverse heuristic.
	MethodPPInv = "partial pseudoinverse"
	// MethodComb is the legacy bounded combinatorial search.
	MethodComb = "combinatorial"
)

// BalanceStrategy selects how [balancer.Auto] builds the coefficient vector.
type BalanceStrategy int

const (
	// StrategyExact is the default: exact nullity-one solving for nullity 1
	// and the exact rational Phase-I simplex for higher nullity. It is sound
	// and complete for every rational reaction matrix.
	StrategyExact BalanceStrategy = iota
	// StrategyLegacy tries the historical heuristics (exact Inv, then GPinv,
	// then PPinv) before falling back to the exact method. It exists for
	// compatibility and comparison; it can return larger valid vectors and is
	// not needed for completeness.
	StrategyLegacy
)

// String returns the strategy name.
func (s BalanceStrategy) String() string {
	switch s {
	case StrategyExact:
		return "exact"
	case StrategyLegacy:
		return "legacy"
	default:
		return "unknown"
	}
}
