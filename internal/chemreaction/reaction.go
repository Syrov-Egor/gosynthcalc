package chemreaction

import (
	"fmt"
	"strings"

	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
	"gonum.org/v1/gonum/mat"
)

//TODO! test case "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2++++"

type ChemicalReaction struct {
	reaction       string
	reacOpts       ReacOptions
	decomposer     *reactionDecomposer
	chemFormulas   *[]chemformula.ChemicalFormula
	parsedFormulas *[][]chemformula.Atom
	molarMasses    *[]float64
	matrix         *mat.Dense
	balancer       *Balancer
	coefs          *MethodResult
	normCoefs      *[]float64
}

type Mode int

const (
	force Mode = iota
	check
	balance
)

type ReacOptions struct {
	mode       Mode
	target     int
	targerMass float64
	intify     bool
	precision  uint
	tolerance  float64
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
			mode:       balance,
			target:     0,
			targerMass: 1.0,
			intify:     true,
			precision:  8,
			tolerance:  1e-8,
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
	if r.reacOpts.target <= high && r.reacOpts.target >= low {
		return r.reacOpts.target - low, nil
	}
	return -1, fmt.Errorf(
		"The target integer %d should be in range %d : %d",
		r.reacOpts.target,
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

func (r *ChemicalReaction) Balancer() (*Balancer, error) {
	mat, err := r.Matrix()
	if err != nil {
		return nil, err
	}
	if r.balancer == nil {
		bal := NewBalancer(
			mat,
			r.decomposer.separatorPos,
			r.reacOpts.intify,
			r.reacOpts.precision,
			r.reacOpts.tolerance)
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
		coeffs := Coeffs{
			mode:               r.reacOpts.mode,
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

		norm = utils.RoundFloatS(norm, r.reacOpts.precision)
		r.normCoefs = &norm
	}
	return *r.normCoefs, nil
}

func (r *ChemicalReaction) isBalanced() bool {
	coefs, _ := r.Coefficients()
	bal, _ := r.Balancer()
	return isReactionBalanced(
		bal.bAlgos.ReactantMatrix,
		bal.bAlgos.ProductMatrix,
		coefs.Result,
		r.reacOpts.tolerance,
	)
}
