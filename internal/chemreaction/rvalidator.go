package chemreaction

import (
	"fmt"
	"strings"
)

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

type reactionValidator struct {
	reaction string
}

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
