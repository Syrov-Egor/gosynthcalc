package chemreaction

import (
	"fmt"
	"strings"
)

// allowedASCII holds the ASCII characters permitted in a reaction string:
// letters, digits, the formula punctuation "., ( { [ ) } ] *" and the
// reaction symbols "=+". allowedExtra holds the few non-ASCII characters
// that are also accepted: the adduct symbols "·" (U+00B7) and "•" (U+2022).
var (
	allowedASCII [128]bool
	allowedExtra = map[rune]bool{
		'·': true, // U+00B7
		'•': true, // U+2022
	}
)

func init() {
	for c := 'a'; c <= 'z'; c++ {
		allowedASCII[c] = true
	}
	for c := 'A'; c <= 'Z'; c++ {
		allowedASCII[c] = true
	}
	for c := '0'; c <= '9'; c++ {
		allowedASCII[c] = true
	}
	for _, c := range ".,({[)}]*=+" {
		allowedASCII[c] = true
	}
}

// sanitize normalizes a reaction string before validation: whitespace is
// dropped and every recognized reactants–products separator (see
// reactionSymbols) is rewritten to the canonical "=" so that the decomposer
// only has to split on a single character.
func sanitize(reaction string) string {
	var res strings.Builder
	res.Grow(len(reaction))
	for _, r := range reaction {
		switch r {
		case ' ':
		default:
			res.WriteRune(r)
		}
	}
	sanitized := res.String()
	for _, sep := range reactionSymbols.reactionSeparators {
		sanitized = strings.ReplaceAll(sanitized, sep, DefaultReactionSeparator)
	}
	return sanitized
}

// reactionValidator checks a raw reaction string before it is decomposed.
type reactionValidator struct {
	reaction string
}

// validate runs the checks in order: the string must not be empty, may
// contain only allowed characters, must carry a reactants–products separator
// (checked while decomposing) and must separate compounds with at least one
// "+". On success it returns the [reactionDecomposer] built from the
// sanitized string.
func (v reactionValidator) validate() (*reactionDecomposer, error) {

	if v.reaction == "" {
		return nil, fmt.Errorf("empty reaction string")
	}

	invalidCharacters := make([]rune, 0)
	for _, r := range v.reaction {
		if r < 128 {
			if !allowedASCII[r] {
				invalidCharacters = append(invalidCharacters, r)
			}
		} else if !allowedExtra[r] {
			invalidCharacters = append(invalidCharacters, r)
		}
	}

	if len(invalidCharacters) > 0 {
		return nil, fmt.Errorf("there are invalid character(s) %s in the reaction '%s'", string(invalidCharacters), v.reaction)
	}

	decomp, err := newReactionDecomposer(v.reaction)
	if err != nil {
		return nil, err
	}

	if !strings.Contains(v.reaction, reactionSymbols.reactantSeparator) {
		return nil, fmt.Errorf("no separators between compounds: %s in the reaction '%s'", reactionSymbols.reactantSeparator, v.reaction)
	}

	return decomp, nil
}
