package chemreaction

import (
	"context"
	"fmt"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

// Auto balances the reaction with the configured [BalanceStrategy] using a
// background context. It is equivalent to AutoContext(context.Background()).
func (b *balancer) Auto() (MethodResult, error) {
	return b.AutoContext(context.Background())
}

// AutoContext balances the reaction exactly and reports the method that
// produced the result.
//
// The default [StrategyExact] decision tree is:
//
//   - nullity 0: no non-zero null vector exists, so the reaction is
//     infeasible ([ErrNoPositiveSolution]);
//   - nullity 1: the unique null direction is checked for zero and mixed-sign
//     entries and returned as primitive integers ([MethodInverse]);
//   - nullity 2 or more: the exact rational Phase-I simplex finds a positive
//     vector or proves infeasibility ([MethodExactLP]).
//
// With [StrategyLegacy], the historical Inv, GPinv and PPinv heuristics are
// tried first, then the exact method is used as the completeness fallback.
// Every legacy candidate is validated as a float with a scaled tolerance.
func (b *balancer) AutoContext(ctx context.Context) (MethodResult, error) {
	if ctx == nil {
		ctx = context.Background()
	}
	if err := ctx.Err(); err != nil {
		return MethodResult{}, err
	}

	a, err := b.analysis()
	if err != nil {
		return MethodResult{}, fmt.Errorf("analyze reaction matrix: %w", err)
	}

	if b.strategy == StrategyLegacy {
		if result, ok := b.tryLegacy(ctx); ok {
			return result, nil
		}
		if err := ctx.Err(); err != nil {
			return MethodResult{}, err
		}
	}

	switch {
	case a.nullity == 0:
		return MethodResult{}, fmt.Errorf("can't balance reaction: %w", ErrNoPositiveSolution)
	case a.nullity == 1:
		coefs, err := b.solveNullityOne(a)
		if err != nil {
			return MethodResult{}, fmt.Errorf("inverse: %w", err)
		}
		return MethodResult{Method: MethodInverse, Result: coefs}, nil
	default:
		coefs, err := b.solveSimplex(ctx, a)
		if err != nil {
			return MethodResult{}, fmt.Errorf("exact LP: %w", err)
		}
		return MethodResult{Method: MethodExactLP, Result: coefs}, nil
	}
}

// tryLegacy runs the historical heuristics in order and returns the first
// validated result. A failed heuristic does not prove infeasibility, so the
// caller still falls back to the exact method.
func (b *balancer) tryLegacy(ctx context.Context) (MethodResult, bool) {
	if err := ctx.Err(); err != nil {
		return MethodResult{}, false
	}
	if coefs, err := b.Inv(); err == nil {
		return MethodResult{Method: MethodInverse, Result: coefs}, true
	}
	if err := ctx.Err(); err != nil {
		return MethodResult{}, false
	}
	if coefs, err := b.GPinv(); err == nil {
		return MethodResult{Method: MethodGPInv, Result: coefs}, true
	}
	if err := ctx.Err(); err != nil {
		return MethodResult{}, false
	}
	if coefs, err := b.PPinv(); err == nil {
		return MethodResult{Method: MethodPPInv, Result: coefs}, true
	}
	return MethodResult{}, false
}

// calculateExactByMethod runs an exact method and returns coefficients that
// were already validated exactly. No rounding or intification is applied.
func (b *balancer) calculateExactByMethod(ctx context.Context, method string) ([]float64, error) {
	errm := fmt.Errorf("can't balance reaction by %s method", method)

	switch method {
	case MethodInverse:
		a, err := b.analysis()
		if err != nil {
			return nil, fmt.Errorf("%w: %w", errm, err)
		}
		coefs, err := b.solveNullityOne(a)
		if err != nil {
			return nil, fmt.Errorf("%w: %w", errm, err)
		}
		return coefs, nil
	case MethodExactLP:
		a, err := b.analysis()
		if err != nil {
			return nil, fmt.Errorf("%w: %w", errm, err)
		}
		coefs, err := b.solveSimplex(ctx, a)
		if err != nil {
			return nil, fmt.Errorf("%w: %w", errm, err)
		}
		return coefs, nil
	default:
		return nil, fmt.Errorf("no exact method %s", method)
	}
}

// calculateLegacyByMethod runs a legacy float heuristic and applies the
// legacy post-processing: rounding to precision+2 decimals, scaled float
// validation and optional intification. maxCoef is read only for
// [MethodComb].
func (b *balancer) calculateLegacyByMethod(ctx context.Context, method string, maxCoef ...uint) ([]float64, error) {
	errm := fmt.Errorf("can't balance reaction by %s method", method)

	var coefficients []float64
	var err error
	switch method {
	case MethodGPInv:
		coefficients, err = b.bAlgos.gPInvAlgorithm()
		if err != nil {
			return nil, fmt.Errorf("%w: %w", errm, err)
		}
	case MethodPPInv:
		coefficients, err = b.bAlgos.pPInvAlgorithm()
		if err != nil {
			return nil, fmt.Errorf("%w: %w", errm, err)
		}
	case MethodComb:
		if len(maxCoef) == 0 {
			return nil, fmt.Errorf("%w: missing maxCoef", errm)
		}
		coefficients = b.bAlgos.combinatorial(ctx, maxCoef[0])
		if coefficients == nil {
			return nil, fmt.Errorf("%w: no combination balances the reaction", errm)
		}
	default:
		return nil, fmt.Errorf("no legacy method %s", method)
	}

	coefficients = utils.RoundFloatS(coefficients, b.precision+2)
	if !b.validFloatResult(coefficients) {
		return nil, fmt.Errorf("wrong coefficients")
	}

	if b.intify {
		intified := b.intifyCoefs(coefficients, b.maxDenom)
		if b.validFloatResult(intified) {
			coefficients = intified
		}
	}
	return coefficients, nil
}
