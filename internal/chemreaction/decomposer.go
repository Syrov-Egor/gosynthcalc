package chemreaction

import (
	"fmt"
	"regexp"
	"strconv"
	"strings"
	"unicode"
)

type regexes struct {
	allowedSymbols     *regexp.Regexp
	reactionSeparators []string
	reactantSeparator  string
}

var reactionRegexes regexes = regexes{
	allowedSymbols: regexp.MustCompile(`[^a-zA-Z0-9.({[)}\]*·•=<\->→⇄+]`),
	reactionSeparators: []string{
		"==",
		"=",
		"<->",
		"->",
		"<>",
		">",
		"→",
		"⇄"},
	reactantSeparator: "+",
}

type compound struct {
	coef    float64
	formula string
}

type rnCoef struct {
	i    int
	coef []rune
}

type reactionDecomposer struct {
	separator    string
	separatorPos int
	initCoefs    []float64
	compounds    []string
	reactants    []string
	products     []string
}

func newReactionDecomposer(reaction string) (*reactionDecomposer, error) {
	if reaction == "" {
		return nil, fmt.Errorf("empty reaction string")
	}

	separator := extractSeparator(reaction)
	initReactants := strings.Split(strings.Split(reaction, separator)[0], reactionRegexes.reactantSeparator)
	initProducts := strings.Split(strings.Split(reaction, separator)[1], reactionRegexes.reactantSeparator)
	splitted := []compound{}
	for i, form := range append(initReactants, initProducts...) {
		if len(form) == 0 {
			return nil, fmt.Errorf("compound %d is empty, maybe there are two adjacent +?", i+1)
		}
		spltCompound, err := splitCoefFromFormula(form)
		if err != nil {
			return nil, err
		}
		splitted = append(splitted, spltCompound)
	}
	initCoefs := make([]float64, len(splitted))
	compounds := make([]string, len(splitted))
	for i, comp := range splitted {
		initCoefs[i] = comp.coef
		compounds[i] = comp.formula
	}

	separatorPos := len(initReactants)

	return &reactionDecomposer{
		separator:    separator,
		separatorPos: separatorPos,
		initCoefs:    initCoefs,
		compounds:    compounds,
		reactants:    compounds[:separatorPos],
		products:     compounds[separatorPos:],
	}, nil
}

func extractSeparator(reaction string) string {
	for _, sep := range reactionRegexes.reactionSeparators {
		if strings.Contains(reaction, sep) {
			splitted := strings.Split(reaction, sep)
			if splitted[0] != "" && splitted[1] != "" {
				return sep
			}
		}
	}
	return ""
}

func splitCoefFromFormula(formula string) (compound, error) {
	if !unicode.IsDigit(rune(formula[0])) {
		return compound{coef: 1.0, formula: formula}, nil
	} else {
		coef := rnCoef{0, []rune{}}
		for i, symbol := range formula {
			if unicode.IsDigit(symbol) || symbol == '.' {
				coef.i = i
				coef.coef = append(coef.coef, symbol)
			} else {
				break
			}
		}
		coefFl, err := strconv.ParseFloat(string(coef.coef), 64)
		return compound{coef: coefFl, formula: formula[coef.i+1:]}, err
	}
}
