package chemreaction

import (
	"fmt"
	"strings"
)

type reactionValidator struct {
	reaction string
}

func (v reactionValidator) emptyReaction() bool {
	return v.reaction == ""
}

func (v reactionValidator) invalidCharacters() []string {
	return reactionRegexes.allowedSymbols.FindAllString(v.reaction, -1)
}

func (v reactionValidator) noRPSeparator(decomp reactionDecomposer) bool {
	return decomp.separator == ""
}

func (v reactionValidator) noReacSeparator() bool {
	return strings.Index(v.reaction, reactionRegexes.reactantSeparator) == -1
}

func (v reactionValidator) validate() (*reactionDecomposer, error) {
	var err error
	decomp, err := NewReactionDecomposer(v.reaction)
	if err != nil {
		return nil, err
	}

	switch {
	case v.emptyReaction():
		err = fmt.Errorf("Empty reaction string")
	case len(v.invalidCharacters()) > 0:
		err = fmt.Errorf("There are invalid character(s) %s in the reaction '%s'",
			v.invalidCharacters(), v.reaction)
	case v.noRPSeparator(*decomp):
		err = fmt.Errorf("No separator between reactants and products: %s in the reaction '%s",
			reactionRegexes.reactionSeparators, v.reaction)
	case v.noReacSeparator():
		err = fmt.Errorf("No separators between compounds: %s in the reaction '%s",
			reactionRegexes.reactantSeparator, v.reaction)
	}

	return decomp, err
}
