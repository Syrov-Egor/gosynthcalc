package chemreaction

import (
	"fmt"
	"strconv"
	"strings"
	"unicode"
)

type symbols struct {
	reactionSeparators []string
	reactantSeparator  string
}

var reactionSymbols symbols = symbols{
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
	separatorPos int
	initCoefs    []float64
	compounds    []string
	reactants    []string
	products     []string
}

func newReactionDecomposer(reaction string) (*reactionDecomposer, error) {
	if len(reaction) < 2 {
		return nil, fmt.Errorf("empty or invalid reaction string")
	}

	reactantsPart, productsPart, found := strings.Cut(reaction, DefaultReactionSeparator)
	if !found || reactantsPart == "" || productsPart == "" {
		return nil, fmt.Errorf("no separator between reactants and products: %s in the reaction '%s'",
			reactionSymbols.reactionSeparators, reaction)
	}
	initReactants := strings.Split(reactantsPart, reactionSymbols.reactantSeparator)
	initProducts := strings.Split(productsPart, reactionSymbols.reactantSeparator)
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
		separatorPos: separatorPos,
		initCoefs:    initCoefs,
		compounds:    compounds,
		reactants:    compounds[:separatorPos],
		products:     compounds[separatorPos:],
	}, nil
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
