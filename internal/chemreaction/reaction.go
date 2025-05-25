package chemreaction

import (
	"strings"

	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"gonum.org/v1/gonum/mat"
)

type ChemicalReaction struct {
	reaction       string
	reacOpts       *ReacOptions
	decomposer     *reactionDecomposer
	chemFormulas   *[]chemformula.ChemicalFormula
	parsedFormulas *[][]chemformula.Atom
	matrix         *mat.Dense
	balancer       *Balancer
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
		reacOpts:   &reacOpt,
	}, nil
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

func (r *ChemicalReaction) Matrix() *mat.Dense {
	if r.matrix == nil {
		parsed, _ := r.ParsedFormulas()
		matrix := createReacMatrix(parsed)
		r.matrix = matrix
	}
	return r.matrix
}

func (r *ChemicalReaction) Balancer() *Balancer {
	if r.balancer == nil {
		bal := NewBalancer(r.Matrix(), r.decomposer.separatorPos, r.reacOpts.intify, r.reacOpts.precision)
		r.balancer = bal
	}
	return r.balancer
}
