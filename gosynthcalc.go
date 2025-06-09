package gosynthcalc

import (
	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"github.com/Syrov-Egor/gosynthcalc/internal/chemreaction"
)

type ChemicalFormula = chemformula.ChemicalFormula

type ChemicalReaction = chemreaction.ChemicalReaction

type ReactionOptions = chemreaction.ReacOptions

type ReactionMode = chemreaction.Mode

const (
	Force   ReactionMode = chemreaction.Force
	Check   ReactionMode = chemreaction.Check
	Balance ReactionMode = chemreaction.Balance
)

func NewChemicalFormula(formula string, precision ...uint) (*ChemicalFormula, error) {
	return chemformula.NewChemicalFormula(formula, precision...)
}

func NewChemicalReaction(reaction string, options ...ReactionOptions) (*ChemicalReaction, error) {
	return chemreaction.NewChemicalReaction(reaction, options...)
}
