package chemreaction

import (
	"bytes"
	"fmt"
	"strconv"
	"strings"
	"text/tabwriter"

	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
	"gonum.org/v1/gonum/mat"
)

type ChemicalReaction struct {
	reaction       string
	reacOpts       ReacOptions
	decomposer     *reactionDecomposer
	chemFormulas   *[]chemformula.ChemicalFormula
	parsedFormulas *[][]chemformula.Atom
	molarMasses    *[]float64
	matrix         *mat.Dense
	balancer       *balancer
	coefs          *MethodResult
	normCoefs      *[]float64
	finalReac      *string
	finalReacNorm  *string
	masses         *[]float64
}

type Mode int

const (
	Force Mode = iota
	Check
	Balance
)

func (m Mode) String() string {
	return [...]string{"Force", "Check", "Balance"}[m]
}

type ReacOptions struct {
	Rmode      Mode
	Target     int
	TargerMass float64
	Intify     bool
	Precision  uint
	Tolerance  float64
}

func NewChemicalReaction(reaction string, options ...ReacOptions) (*ChemicalReaction, error) {
	newReaction := strings.Replace(reaction, " ", "", -1)
	validator := reactionValidator{reaction: newReaction}
	decomp, err := validator.validate()
	if err != nil {
		return nil, err
	}

	var reacOpt ReacOptions
	if options == nil {
		reacOpt = ReacOptions{
			Rmode:      Balance,
			Target:     0,
			TargerMass: 1.0,
			Intify:     true,
			Precision:  8,
			Tolerance:  1e-8,
		}
	} else {
		reacOpt = options[0]
	}

	return &ChemicalReaction{
		reaction:   newReaction,
		decomposer: decomp,
		reacOpts:   reacOpt,
	}, nil
}

func (r *ChemicalReaction) calculatedTarget() (int, error) {
	high := len(r.decomposer.products) - 1
	low := -len(r.decomposer.reactants)
	if r.reacOpts.Target <= high && r.reacOpts.Target >= low {
		return r.reacOpts.Target - low, nil
	}
	return -1, fmt.Errorf(
		"The target integer %d should be in range %d : %d",
		r.reacOpts.Target,
		low,
		high,
	)
}

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

func (r *ChemicalReaction) Balancer() (*balancer, error) {
	mat, err := r.Matrix()
	if err != nil {
		return nil, err
	}
	if r.balancer == nil {
		bal := newBalancer(
			mat,
			r.decomposer.separatorPos,
			r.reacOpts.Intify,
			r.reacOpts.Precision,
			r.reacOpts.Tolerance)
		r.balancer = bal
	}
	return r.balancer, nil
}

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

func (r *ChemicalReaction) SetCoefficients(coefs []float64) error {
	if len(coefs) != len(r.decomposer.compounds) {
		return fmt.Errorf("Lenght of coefficient slice should be %d, got %d", len(r.decomposer.compounds), len(coefs))
	}
	for i, coef := range coefs {
		if coef <= 0 {
			return fmt.Errorf("Input coefficient %f at position %d is <= 0", coef, i)
		}
	}

	r.coefs = &MethodResult{Method: "User", Result: coefs}

	return nil
}

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

func (r *ChemicalReaction) IsBalanced() bool {
	coefs, _ := r.Coefficients()
	bal, _ := r.Balancer()
	return isReactionBalanced(
		bal.bAlgos.ReactantMatrix,
		bal.bAlgos.ProductMatrix,
		coefs.Result,
		r.reacOpts.Tolerance,
	)
}

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
		reactionRegexes.reactantSeparator,
		r.decomposer.separator,
		r.decomposer.separatorPos,
	)

	return replaced
}

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
		nu := r.reacOpts.TargerMass / molars[target]
		masses := make([]float64, len(molars))
		for i, molar := range molars {
			masses[i] = utils.RoundFloat(molar*nu*normCoefs[i], r.reacOpts.Precision)
		}

		r.masses = &masses
	}
	return *r.masses, nil
}

func (r *ChemicalReaction) Output(printPrecision ...uint) (crOutput, error) {
	var pPrecision uint
	if printPrecision == nil {
		pPrecision = 4
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
		MolarMasses:       utils.RoundFloatS(mMasses, pPrecision),
		Target:            r.decomposer.compounds[target],
		Masses:            utils.RoundFloatS(mass, pPrecision),
	}
	return crO, nil
}

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
}

func (o crOutput) String() string {
	out := fmt.Sprintln("initial reaction:", o.Reaction) +
		fmt.Sprint("reaction matrix:\n", o.Matrix, "\n") +
		fmt.Sprintln("mode:", o.Mode) +
		fmt.Sprintln("formulas:", o.Formulas) +
		fmt.Sprintln("coefficients:", o.Coefficients) +
		fmt.Sprintln("coefficients normalized:", o.NormCoefficients) +
		fmt.Sprintln("algorithm:", o.Algorithm) +
		fmt.Sprintln("is balanced:", o.IsBalanced) +
		fmt.Sprintln("final reaction:", o.FinalReaction) +
		fmt.Sprintln("final reaction normalized:", o.FinalReactionNorm) +
		fmt.Sprintln("molar masses:", o.MolarMasses) +
		fmt.Sprintln("target:", o.Target) +
		fmt.Sprintln("masses:", o.Masses)

	var buf bytes.Buffer
	w := tabwriter.NewWriter(&buf, 0, 0, 2, ' ', 0)

	for i, comp := range o.Formulas {
		fmt.Fprintf(w, "%s\tM = %v\tg/mol\tm = %v\tg\n",
			comp, o.MolarMasses[i], o.Masses[i])
	}

	w.Flush()
	return out + strings.TrimSuffix(buf.String(), "\n")
}
