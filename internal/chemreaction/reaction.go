package chemreaction

import (
	"strings"

	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"gonum.org/v1/gonum/mat"
)

type ChemicalReaction struct {
	reaction       string
	decomposer     *reactionDecomposer
	chemFormulas   *[]chemformula.ChemicalFormula
	parsedFormulas *[][]chemformula.Atom
	matrix         *mat.Dense
}

type ReacOptions struct {
}

func NewChemicalReaction(reaction string, options ...ReacOptions) (*ChemicalReaction, error) {
	newReaction := strings.Replace(reaction, " ", "", -1)
	validator := reactionValidator{reaction: newReaction}
	decomp, err := validator.validate()
	if err != nil {
		return nil, err
	}

	return &ChemicalReaction{
		reaction:   newReaction,
		decomposer: decomp,
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
