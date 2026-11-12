// Package chemreaction parses, balances and evaluates chemical reactions.
//
// The entry point is [ChemicalReaction]: a reaction
// string such as "KIO3+KI+H2SO4=I2+K2SO4+H2O" is decomposed into reactants
// and products, assembled into a reaction matrix (rows are elements, columns
// are compounds) and then balanced, either by checking or reusing the
// coefficients written in the string or by the exact or legacy algorithms
// exposed through [balancer]. From the coefficients the package derives the
// normalized coefficients, the final reaction strings and, finally, the
// masses of every compound for a chosen target and target mass.
//
// # Modes
//
// Three coefficient calculation modes are available (see [Mode]):
//
//  1. Force uses the coefficients entered in the reaction string and
//     calculates the masses whether the reaction is balanced or not.
//
//  2. Check is the same as Force, but verifies that the entered
//     coefficients balance the reaction.
//
//  3. Balance (the default) calculates the coefficients automatically.
//
// # Coefficient algorithms
//
// The default [StrategyExact] uses exact rational arithmetic and classifies
// the signed reaction matrix A = [R | -P] by its nullity:
//
//   - nullity 0 means A x = 0 has only x = 0; the reaction is infeasible and
//     the error wraps [ErrNoPositiveSolution];
//   - nullity 1 gives a line span{v}; v is accepted only if it has no zero
//     entries and all entries have the same sign, then scaled to primitive
//     integers (method [MethodInverse]);
//   - nullity 2 or more is solved by the exact rational Phase-I simplex
//     (method [MethodExactLP]), which either finds a positive vector or
//     proves infeasibility.
//
// Nullity alone is not an existence test: a nullity-1 direction can contain a
// zero or mixed signs, and a higher-nullity space can still contain no
// positive vector. Both cases return [ErrNoPositiveSolution].
//
// [StrategyLegacy] tries the historical Inv, GPinv and PPinv heuristics first
// and then falls back to the exact method. It exists for compatibility and
// comparison; its results can be larger and it is not needed for completeness.
// A failed heuristic returns [ErrMethodNotPositive]; it never claims that the
// reaction is infeasible. Intify, Precision and Tolerance affect only legacy
// float results.
//
// Exact coefficients are primitive positive integers and are never rounded or
// intified. [ChemicalReaction.CoefficientsExact] returns them without passing
// through float64, so values above 2^53 are preserved; the float API returns
// an error wrapping [ErrCoefficientOverflow] instead of a corrupted value.
// [ChemicalReaction.NormCoefficients] still divides by the target coefficient
// and rounds to Precision for display. The simplex also returns
// [ErrSimplexIterationLimit] on its safety limit and propagates context
// cancellation unchanged.
//
// # Accepted input
//
// A reaction string may use any of the separators "==", "=", "<->", "->",
// "<>", ">", "→" or "⇄" between reactants and products (all normalized to
// "="), "+" between compounds, and leading numeric coefficients such as
// "3H2SO4". Whitespace is ignored. The accepted characters are letters,
// digits, the formula punctuation ". ( { [ ) } ] *", the reaction separators
// and "+"; anything else makes validation fail. For anything but Force mode
// the atom sets of the two sides of the equation must match, so that a
// solution of the equation exists at all.
package chemreaction

import (
	"bytes"
	"context"
	"fmt"
	"math/big"
	"slices"
	"strconv"
	"strings"
	"text/tabwriter"

	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
	"gonum.org/v1/gonum/mat"
)

const (
	// DefaultMode is the coefficients calculation mode used when no options
	// are supplied to [NewChemicalReaction]: [Balance].
	DefaultMode Mode = Balance
	// DefaultTarget is the default target compound: 0, the first product.
	DefaultTarget int = 0
	// DefaultTargetMass is the default mass of the target compound in grams.
	DefaultTargetMass float64 = 1.0
	// DefaultIntify is the default value of [ReacOptions.Intify]: integer
	// coefficients are preferred.
	DefaultIntify bool = true
	// DefaultPrecision is the default rounding precision.
	DefaultPrecision uint = 8
	// DefaultPrintPrecision is the default number of digits printed by
	// [ChemicalReaction.Output].
	DefaultPrintPrecision uint = 4
	// DefaultTolerance is the default float tolerance for legacy comparisons,
	// used when [ReacOptions.Tolerance] is not set. Balance checks scale it by
	// the coefficient magnitude.
	DefaultTolerance float64 = 1e-8
	// DefaultReactionSeparator is the canonical separator between reactants
	// and products ("="). Every separator listed in reactionSymbols is
	// normalized to it during validation, and generated final reaction
	// strings always use it.
	DefaultReactionSeparator string = "="
)

// ChemicalReaction represents a chemical reaction string and computes
// everything that can be derived from it: the chemical formulas of its
// compounds, their parsed forms and molar masses, the reaction matrix, the
// coefficients (taken from the string, checked, or balanced automatically),
// their normalization onto the target compound, the final reaction strings
// and the masses of all compounds. Construct it with
// [NewChemicalReaction].
//
// Every property is calculated lazily on first use and cached afterwards;
// [ChemicalReaction.SetCoefficients] invalidates the caches that depend on
// the coefficients. Slices returned by the accessors are the cached values
// themselves: callers must not modify them.
type ChemicalReaction struct {
	reaction            string
	reacOpts            ReacOptions
	decomposer          *reactionDecomposer
	chemFormulas        *[]chemformula.ChemicalFormula
	parsedFormulas      *[][]chemformula.Atom
	exactParsedFormulas *[][]chemformula.ExactAtom
	molarMasses         *[]float64
	matrix              *mat.Dense
	balancer            *balancer
	coefs               *MethodResult
	normCoefs           *[]float64
	finalReac           *string
	finalReacNorm       *string
	masses              *[]float64
}

// Mode selects how the coefficients of a reaction are obtained.
type Mode int

const (
	// Force uses the coefficients written in the reaction string exactly as
	// they are, without any balance checking: masses are calculated whether
	// the reaction is balanced or not. This mode does not imply any
	// automatic checking of calculation results.
	Force Mode = iota
	// Check is like [Force], but the coefficients must balance the reaction;
	// otherwise the calculation fails with "reaction is not balanced".
	Check
	// Balance tries to calculate the coefficients from the reaction string
	// automatically with [balancer.Auto]: exact nullity-one solving or the
	// exact rational Phase-I simplex by default, or the historical heuristics
	// first when [ReacOptions.Strategy] is [StrategyLegacy]. It is the default
	// mode.
	Balance
)

// String returns the name of the mode ("Force", "Check" or "Balance").
// Calling it with a value outside the defined modes panics.
func (m Mode) String() string {
	return [...]string{"Force", "Check", "Balance"}[m]
}

// ReacOptions configures a [ChemicalReaction]. When [NewChemicalReaction] is
// called without options, the Default* constants of this package are used;
// when options are passed, they are taken as-is, so every field should be
// filled in.
type ReacOptions struct {
	// Rmode is the coefficients calculation mode: [Force], [Check] or
	// [Balance] (the default).
	Rmode Mode
	// Target is the index of the target compound, i.e. the substance whose
	// mass is known in advance: 0 selects the first product (the usual
	// choice), positive values walk forward through the products, and
	// negative values select reactants counted from the right, -1 being the
	// last reactant (the usual limiting reagent) and -number of reactants
	// the first one. Out-of-range values fail in [ChemicalReaction.Masses]
	// and friends with "the target integer %d should be in range %d : %d".
	Target int
	// TargetMass is the desired mass of the target compound in grams.
	TargetMass float64
	// Intify tells the legacy float path to convert fractional coefficients
	// into the smallest integers with the same ratios; when it cannot be done
	// safely the floats are kept. The exact path already returns primitive
	// integers and ignores Intify. Coefficients are never intified in [Force]
	// and [Check] modes, which use the numbers from the reaction string.
	Intify bool
	// Precision is the number of decimal places used when rounding calculated
	// values (molar masses, percentages, normalized coefficients, masses) and
	// when rounding legacy float coefficients. Exact coefficients are never
	// rounded.
	Precision uint
	// Tolerance is the absolute float tolerance used for legacy candidate
	// validation and [Check] mode; comparisons scale it by the coefficient
	// magnitude. It does not affect the exact path or [ChemicalReaction.IsBalanced]
	// for integer coefficients.
	Tolerance float64
	// Strategy selects how [Mode] Balance computes coefficients. The default
	// zero value, [StrategyExact], uses the exact nullity-one solver and the
	// exact rational Phase-I simplex. [StrategyLegacy] tries the historical
	// heuristics first and falls back to the exact method.
	Strategy BalanceStrategy
}

// NewChemicalReaction validates and decomposes reaction and returns a
// [ChemicalReaction] configured with options. When no options are given, the
// Default* constants of this package are used; otherwise options[0] is taken
// as-is.
//
// Validation rejects an empty string, characters outside the accepted set
// (see the package documentation), a missing reactants–products separator, a
// missing "+" between compounds and empty compounds (two adjacent "+").
//
// Example:
//
//	reaction, err := NewChemicalReaction("H2+O2=H2O")
//	reaction.Coefficients()
func NewChemicalReaction(reaction string, options ...ReacOptions) (*ChemicalReaction, error) {
	sanReaction := sanitize(reaction)
	validator := reactionValidator{reaction: sanReaction}
	decomp, err := validator.validate()
	if err != nil {
		return nil, err
	}

	var reacOpt ReacOptions
	if options == nil {
		reacOpt = ReacOptions{
			Rmode:      DefaultMode,
			Target:     DefaultTarget,
			TargetMass: DefaultTargetMass,
			Intify:     DefaultIntify,
			Precision:  DefaultPrecision,
			Tolerance:  DefaultTolerance,
			Strategy:   StrategyExact,
		}
	} else {
		reacOpt = options[0]
	}

	return &ChemicalReaction{
		reaction:   reaction,
		decomposer: decomp,
		reacOpts:   reacOpt,
	}, nil
}

// calculatedTarget maps [ReacOptions.Target] onto a plain index into the
// flat list of compounds (reactants first, then products) and reports an
// error when the target lies outside the reaction.
func (r *ChemicalReaction) calculatedTarget() (int, error) {
	high := len(r.decomposer.products) - 1
	low := -len(r.decomposer.reactants)
	if r.reacOpts.Target <= high && r.reacOpts.Target >= low {
		return r.reacOpts.Target - low, nil
	}
	return -1, fmt.Errorf(
		"the target integer %d should be in range %d : %d",
		r.reacOpts.Target,
		low,
		high,
	)
}

// ChemFormulas returns a [chemformula.ChemicalFormula] for every compound of
// the reaction, reactants first. An error is returned if some compound
// fails to parse on its own, e.g. because it contains an unknown element
// symbol or unbalanced parentheses that the reaction-level validation does not catch.
func (r *ChemicalReaction) ChemFormulas() ([]chemformula.ChemicalFormula, error) {
	if r.chemFormulas == nil {
		formulas := []chemformula.ChemicalFormula{}
		for _, compound := range r.decomposer.compounds {
			f, err := chemformula.NewChemicalFormula(compound)
			if err != nil {
				return nil, err
			}
			formulas = append(formulas, *f)
		}
		r.chemFormulas = &formulas
	}
	return *r.chemFormulas, nil
}

// ParsedFormulas returns the parsed form of every compound, each in the
// order its elements first appear in the formula.
func (r *ChemicalReaction) ParsedFormulas() ([][]chemformula.Atom, error) {
	if r.parsedFormulas == nil {
		c, err := r.ChemFormulas()
		if err != nil {
			return nil, err
		}
		parsed := [][]chemformula.Atom{}
		for _, compound := range c {
			parsed = append(parsed, compound.ParsedFormula())
		}

		r.parsedFormulas = &parsed
	}

	return *r.parsedFormulas, nil
}

// ExactParsedFormulas returns the exact rational parsed form of every
// compound, each in the order its elements first appear in the formula. It is
// the exact counterpart of [ChemicalReaction.ParsedFormulas] and avoids
// binary floating-point artifacts from nested fractional groups. The returned
// slices are cached; callers must not modify them.
func (r *ChemicalReaction) ExactParsedFormulas() ([][]chemformula.ExactAtom, error) {
	if r.exactParsedFormulas == nil {
		c, err := r.ChemFormulas()
		if err != nil {
			return nil, err
		}
		parsed := [][]chemformula.ExactAtom{}
		for _, compound := range c {
			parsed = append(parsed, compound.ExactParsedFormula())
		}
		r.exactParsedFormulas = &parsed
	}
	return *r.exactParsedFormulas, nil
}

// MolarMasses returns the [chemformula.ChemicalFormula.MolarMass] of every
// compound in g/mol.
func (r *ChemicalReaction) MolarMasses() ([]float64, error) {
	if r.molarMasses == nil {
		masses := make([]float64, len(r.decomposer.compounds))
		formulas, err := r.ChemFormulas()
		if err != nil {
			return nil, err
		}
		for i, formula := range formulas {
			masses[i] = formula.MolarMass()
		}
		r.molarMasses = &masses
	}
	return *r.molarMasses, nil
}

// Matrix returns the reaction matrix: rows are the distinct elements in the
// order of their first appearance in the reaction, columns are the
// compounds. See [createReacMatrix] for the construction and the Blakley
// reference.
func (r *ChemicalReaction) Matrix() (*mat.Dense, error) {
	if r.matrix == nil {
		parsed, err := r.ParsedFormulas()
		if err != nil {
			return nil, err
		}
		matrix := createReacMatrix(parsed)
		r.matrix = matrix
	}
	return r.matrix, nil
}

// Balancer returns the balancer of this reaction: the object that exposes
// the individual balancing algorithms ([balancer.Inv], [balancer.GPinv],
// [balancer.PPinv], [balancer.Comb]), the exact universal solver and
// [balancer.Auto]. It is configured with the reaction matrix, the position of
// the separator, Intify, Precision, Strategy and Tolerance from [ReacOptions].
// Note that the balancer works in every mode: coefficients computed through it
// can be assigned with [ChemicalReaction.SetCoefficients].
func (r *ChemicalReaction) Balancer() (*balancer, error) {
	mat, err := r.Matrix()
	if err != nil {
		return nil, err
	}
	if r.balancer == nil {
		exactParsed, err := r.ExactParsedFormulas()
		if err != nil {
			return nil, err
		}
		exact, err := createExactReacMatrix(exactParsed, r.decomposer.separatorPos)
		if err != nil {
			return nil, err
		}
		bal := newBalancerWithExact(
			mat,
			exact,
			r.decomposer.separatorPos,
			r.reacOpts.Intify,
			r.reacOpts.Precision,
			r.reacOpts.Strategy,
			r.reacOpts.Tolerance)
		r.balancer = bal
	}
	return r.balancer, nil
}

// Coefficients returns the coefficients of the reaction together with the
// name of the method that produced them. Depending on [ReacOptions.Rmode]
// they are taken from the reaction string, checked against the balance condition, or computed
// automatically; the element sets of both sides are verified first in every
// mode except [Force]. The result is cached until [ChemicalReaction.SetCoefficients]
// replaces it.
func (r *ChemicalReaction) Coefficients() (*MethodResult, error) {
	if r.coefs == nil {
		parsed, err := r.ParsedFormulas()
		if err != nil {
			return nil, err
		}
		bal, err := r.Balancer()
		if err != nil {
			return nil, err
		}
		coeffs := coeffs{
			mode:               r.reacOpts.Rmode,
			parsedFormulas:     parsed,
			decomposedReaction: r.decomposer,
			balancer:           bal,
		}
		coefs, err := coeffs.getCoeffs()
		if err != nil {
			return nil, err
		}
		r.coefs = &coefs
	}
	return r.coefs, nil
}

// ExactMethodResult is a primitive integer coefficient list together with the
// name of the method that produced it. It is the full-precision counterpart of
// [MethodResult] and is not limited by the exact float64 range.
type ExactMethodResult struct {
	// Method is the name of the algorithm or source that produced Result.
	Method string
	// Result holds one positive primitive integer per compound.
	Result []*big.Int
}

// CoefficientsExact returns the coefficients as primitive positive integers
// without passing through float64. For the default exact strategy it calls the
// exact solver directly, so coefficients above 2^53 are preserved. For
// user-supplied, Force, Check or legacy coefficients it converts the current
// float values to their exact decimal rationals and scales them to primitive
// integers; the ratios are preserved but the display scale may differ.
//
// An infeasible automatic balance returns an error wrapping
// [ErrNoPositiveSolution].
func (r *ChemicalReaction) CoefficientsExact() (*ExactMethodResult, error) {
	useExactSolver := r.reacOpts.Rmode == Balance &&
		r.reacOpts.Strategy == StrategyExact &&
		(r.coefs == nil || r.coefs.Method != MethodUser)
	if !useExactSolver {
		coefs, err := r.Coefficients()
		if err != nil {
			return nil, err
		}
		return exactResultFromFloats(coefs)
	}

	bal, err := r.Balancer()
	if err != nil {
		return nil, err
	}
	a, err := bal.analysis()
	if err != nil {
		return nil, fmt.Errorf("analyze reaction matrix: %w", err)
	}

	switch {
	case a.nullity == 0:
		return nil, fmt.Errorf("can't balance reaction: %w", ErrNoPositiveSolution)
	case a.nullity == 1:
		z, err := bal.solveNullityOneExact(a)
		if err != nil {
			return nil, fmt.Errorf("inverse: %w", err)
		}
		return &ExactMethodResult{Method: MethodInverse, Result: z}, nil
	default:
		z, err := bal.solveSimplexExact(context.Background(), a)
		if err != nil {
			return nil, fmt.Errorf("exact LP: %w", err)
		}
		return &ExactMethodResult{Method: MethodExactLP, Result: z}, nil
	}
}

// exactResultFromFloats converts a float method result into primitive integer
// ratios.
func exactResultFromFloats(coefs *MethodResult) (*ExactMethodResult, error) {
	if coefs == nil {
		return nil, fmt.Errorf("%w: nil coefficient result", ErrInvalidMatrix)
	}
	x := make([]*big.Rat, len(coefs.Result))
	for i, f := range coefs.Result {
		r, err := ratFromFloat(f)
		if err != nil {
			return nil, fmt.Errorf("coefficient %d: %w", i, err)
		}
		x[i] = r
	}
	z, err := primitiveIntegers(x)
	if err != nil {
		return nil, err
	}
	return &ExactMethodResult{Method: coefs.Method, Result: z}, nil
}

// SetCoefficients replaces the coefficients with a user-supplied list. The
// stored [MethodResult] then carries the method name [MethodUser]. The list must
// have one entry per compound and every entry must be positive. Setting new
// coefficients invalidates the cached normalized coefficients, final reaction
// strings and masses, which are then recomputed from them.
func (r *ChemicalReaction) SetCoefficients(coefs []float64) error {
	if len(coefs) != len(r.decomposer.compounds) {
		return fmt.Errorf("length of coefficient slice should be %d, got %d", len(r.decomposer.compounds), len(coefs))
	}
	for i, coef := range coefs {
		if coef <= 0 {
			return fmt.Errorf("input coefficient %f at position %d is <= 0", coef, i)
		}
	}

	r.coefs = &MethodResult{Method: MethodUser, Result: slices.Clone(coefs)}
	r.normCoefs = nil
	r.finalReac = nil
	r.finalReacNorm = nil
	r.masses = nil

	return nil
}

// NormCoefficients returns the coefficients normalized onto the target
// compound, so that the coefficient of the target itself equals 1, rounded to
// [ReacOptions.Precision].
func (r *ChemicalReaction) NormCoefficients() ([]float64, error) {
	if r.normCoefs == nil {
		coefs, err := r.Coefficients()
		if err != nil {
			return nil, err
		}
		calc, err := r.calculatedTarget()
		if err != nil {
			return nil, err
		}
		targetCompound := coefs.Result[calc]
		norm := make([]float64, len(coefs.Result))
		for i, coef := range coefs.Result {
			norm[i] = coef / targetCompound
		}

		norm = utils.RoundFloatS(norm, r.reacOpts.Precision)
		r.normCoefs = &norm
	}
	return *r.normCoefs, nil
}

// IsBalanced reports whether the current coefficients balance the reaction.
// Integer coefficients are verified exactly against the rational reaction
// matrix; fractional coefficients are checked with a float tolerance scaled
// by the coefficient magnitude. Missing coefficients, a wrong number of
// coefficients or any other failure on the way makes it return false rather
// than an error.
func (r *ChemicalReaction) IsBalanced() bool {
	coefs, err := r.Coefficients()
	if err != nil || coefs == nil {
		return false
	}
	bal, err := r.Balancer()
	if err != nil || bal == nil || bal.bAlgos == nil {
		return false
	}

	reactantMatrix := bal.bAlgos.ReactantMatrix
	productMatrix := bal.bAlgos.ProductMatrix
	if reactantMatrix == nil || productMatrix == nil {
		return false
	}

	_, reactantCols := reactantMatrix.Dims()
	_, productCols := productMatrix.Dims()
	if len(coefs.Result) != reactantCols+productCols {
		return false
	}

	// Integer coefficients can be checked exactly. This is stronger than a
	// scaled float tolerance and is not fooled by a large configured
	// tolerance. Fractional user coefficients fall back to float validation.
	if A, _, analysisErr := bal.exactMatrixAndAnalysis(); analysisErr == nil {
		if z, ok := exactIntegerCoefficients(coefs.Result); ok {
			return validExactResult(A, z)
		}
	}

	return isReactionBalanced(
		reactantMatrix,
		productMatrix,
		coefs.Result,
		scaledTolerance(coefs.Result, r.reacOpts.Tolerance),
	)
}

// generateFinalReaction renders the reaction string from coefs: a
// coefficient equal to 1 is omitted, the compounds are joined with "+" and
// the separatorPos-th "+" — the one between the last reactant and the first
// product — is replaced by [DefaultReactionSeparator], as in
// "2KMnO4+16HCl=2MnCl2+5Cl2+8H2O+2KCl".
func (r *ChemicalReaction) generateFinalReaction(coefs []float64) string {
	final := []string{}
	for i, compound := range r.decomposer.compounds {
		if coefs[i] != 1.0 {
			final = append(final, strconv.FormatFloat(
				coefs[i],
				'f',
				-1,
				64))
		}
		final = append(final, compound)
		final = append(final, "+")
	}
	joined := strings.Join(final[:len(final)-1], "")
	replaced := utils.ReplaceNthOccurrence(
		joined,
		reactionSymbols.reactantSeparator,
		DefaultReactionSeparator,
		r.decomposer.separatorPos,
	)

	return replaced
}

// FinalReaction returns the reaction string with the calculated
// coefficients, e.g. "2KMnO4+16HCl=2MnCl2+5Cl2+8H2O+2KCl".
func (r *ChemicalReaction) FinalReaction() (string, error) {
	if r.finalReac == nil {
		coefs, err := r.Coefficients()
		if err != nil {
			return "", err
		}
		fin := r.generateFinalReaction(coefs.Result)
		r.finalReac = &fin
	}
	return *r.finalReac, nil
}

// FinalReactionNorm returns the reaction string with the normalized
// coefficients, the target therefore carrying no number, e.g.
// "KMnO4+8HCl=MnCl2+2.5Cl2+4H2O+KCl".
func (r *ChemicalReaction) FinalReactionNorm() (string, error) {
	if r.finalReacNorm == nil {
		coefs, err := r.NormCoefficients()
		if err != nil {
			return "", err
		}
		fin := r.generateFinalReaction(coefs)
		r.finalReacNorm = &fin
	}
	return *r.finalReacNorm, nil
}

// Masses returns the mass in grams of every compound needed to obtain
// [ReacOptions.TargetMass] grams of the target compound. The amount of
// substance of the target, nu = target mass / molar mass, is broadcast
// to the other compounds through their normalized coefficients:
// mass = molar mass * nu * normalized coefficient, with each
// value rounded to [ReacOptions.Precision].
func (r *ChemicalReaction) Masses() ([]float64, error) {
	if r.masses == nil {
		molars, err := r.MolarMasses()
		if err != nil {
			return nil, err
		}
		target, err := r.calculatedTarget()
		if err != nil {
			return nil, err
		}
		normCoefs, err := r.NormCoefficients()
		if err != nil {
			return nil, err
		}
		nu := r.reacOpts.TargetMass / molars[target]
		masses := make([]float64, len(molars))
		for i, molar := range molars {
			masses[i] = utils.RoundFloat(molar*nu*normCoefs[i], r.reacOpts.Precision)
		}

		r.masses = &masses
	}
	return *r.masses, nil
}

// Output collects every computed property of the reaction into a [crOutput]
// value in a single pass. The optional printPrecision overrides
// [DefaultPrintPrecision] and controls only how many digits [crOutput.String]
// prints, not the stored values. An error is returned if any property fails
// to compute (bad target, unbalanced reaction, impossible coefficients).
func (r *ChemicalReaction) Output(printPrecision ...uint) (crOutput, error) {
	var pPrecision uint
	if printPrecision == nil {
		pPrecision = DefaultPrintPrecision
	} else {
		pPrecision = printPrecision[0]
	}

	matr, err := r.Matrix()
	if err != nil {
		return crOutput{}, err
	}
	matrix := mat.Formatted(matr)
	coefs, err := r.Coefficients()
	if err != nil {
		return crOutput{}, err
	}
	ncoefs, err := r.NormCoefficients()
	if err != nil {
		return crOutput{}, err
	}
	fReaction, err := r.FinalReaction()
	if err != nil {
		return crOutput{}, err
	}
	nfReaction, err := r.FinalReactionNorm()
	if err != nil {
		return crOutput{}, err
	}
	mMasses, err := r.MolarMasses()
	if err != nil {
		return crOutput{}, err
	}
	target, err := r.calculatedTarget()
	if err != nil {
		return crOutput{}, err
	}
	mass, err := r.Masses()
	if err != nil {
		return crOutput{}, err
	}

	crO := crOutput{
		Reaction:          r.reaction,
		Matrix:            fmt.Sprintf("%v", matrix),
		Mode:              r.reacOpts.Rmode.String(),
		Formulas:          r.decomposer.compounds,
		Coefficients:      coefs.Result,
		NormCoefficients:  ncoefs,
		Algorithm:         coefs.Method,
		IsBalanced:        r.IsBalanced(),
		FinalReaction:     fReaction,
		FinalReactionNorm: nfReaction,
		MolarMasses:       mMasses,
		Target:            r.decomposer.compounds[target],
		Masses:            mass,
		printPrecision:    pPrecision,
	}
	return crO, nil
}

// crOutput is the collected result of a [ChemicalReaction] computation. Its
// [crOutput.String] method renders the same human-readable report, including
// the trailing table of molar masses and masses per compound.
type crOutput struct {
	Reaction          string
	Matrix            string
	Mode              string
	Formulas          []string
	Coefficients      []float64
	NormCoefficients  []float64
	Algorithm         string
	IsBalanced        bool
	FinalReaction     string
	FinalReactionNorm string
	MolarMasses       []float64
	Target            string
	Masses            []float64
	printPrecision    uint
}

// String formats the report line by line: the initial reaction, the reaction
// matrix, the mode, the formulas, the coefficients and their
// normalization, the algorithm, the balance flag, the final reaction
// strings, the molar masses, the target and the masses, followed by a
// tab-aligned "M = ... g/mol, m = ... g" table. Floats are printed with
// printPrecision digits.
func (o crOutput) String() string {
	var out strings.Builder
	fmt.Fprintln(&out, "initial reaction:", o.Reaction)
	fmt.Fprint(&out, "reaction matrix:\n", o.Matrix, "\n")
	fmt.Fprintln(&out, "mode:", o.Mode)
	fmt.Fprintln(&out, "formulas:", o.Formulas)
	fmt.Fprintln(&out, "coefficients:", o.Coefficients)
	fmt.Fprintln(&out, "coefficients normalized:", o.NormCoefficients)
	fmt.Fprintln(&out, "algorithm:", o.Algorithm)
	fmt.Fprintln(&out, "is balanced:", o.IsBalanced)
	fmt.Fprintln(&out, "final reaction:", o.FinalReaction)
	fmt.Fprintln(&out, "final reaction normalized:", o.FinalReactionNorm)
	fmt.Fprintln(&out, "molar masses:", formatFloats(o.MolarMasses, o.printPrecision))
	fmt.Fprintln(&out, "target:", o.Target)
	fmt.Fprintln(&out, "masses:", formatFloats(o.Masses, o.printPrecision))

	var buf bytes.Buffer
	w := tabwriter.NewWriter(&buf, 0, 0, 2, ' ', 0)

	for i, comp := range o.Formulas {
		_, _ = fmt.Fprintf(w, "%s\tM = %.*f\tg/mol\tm = %.*f\tg\n",
			comp, o.printPrecision, o.MolarMasses[i], o.printPrecision, o.Masses[i])
	}
	_ = w.Flush()

	return out.String() + strings.TrimSuffix(buf.String(), "\n")
}

// formatFloats renders a slice of floats as "[1.00000000 8.00000000 ...]"
// with every value printed with precision digits.
func formatFloats(vals []float64, precision uint) string {
	var sb strings.Builder
	sb.WriteByte('[')
	for i, val := range vals {
		if i > 0 {
			sb.WriteByte(' ')
		}
		fmt.Fprintf(&sb, "%.*f", precision, val)
	}
	sb.WriteByte(']')
	return sb.String()
}
