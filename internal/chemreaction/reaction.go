package chemreaction

import "fmt"

type ChemicalReaction struct {
	reaction   string
	decomposer *reactionDecomposer
}

type ReacOptions struct {
}

func NewChemicalReaction(reaction string, options ...ReacOptions) (*ChemicalReaction, error) {
	decomposer, err := NewReactionDecomposer(reaction)
	fmt.Println(decomposer)
	return nil, err
}
