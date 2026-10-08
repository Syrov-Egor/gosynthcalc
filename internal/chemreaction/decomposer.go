package chemreaction

import (
	"fmt"
	"strconv"
	"strings"
	"unicode"
)

// symbols collects the textual conventions of a reaction string: the
// separators allowed between reactants and products, and the only separator
// allowed between compounds.
type symbols struct {
	reactionSeparators []string
	reactantSeparator  string
}

// reactionSymbols lists every separator accepted in a reaction string
// ("==", "=", "<->", "->", "<>", ">", "→", "⇄") plus the compound separator
// "+". All reaction separators are normalized to "=" by [sanitize] before
// decomposition.
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

// compound pairs a formula string with the numeric coefficient that was
// written in front of it in the reaction string (1.0 when there was none).
type compound struct {
	coef    float64
	formula string
}

// rnCoef is the scratch state of [splitCoefFromFormula]: i is the index of
// the last digit or decimal point scanned and coef collects those runes.
type rnCoef struct {
	i    int
	coef []rune
}

// reactionDecomposer splits a reaction string into reactants and products and
// strips the numeric coefficients off the formulas.
//
// Fields:
//
//   - separatorPos is the number of reactants, i.e. the column at which the
//     products begin in the reaction matrix;
//   - initCoefs are the coefficients found in the reaction string (1.0 where
//     none were written);
//   - compounds is every formula without its coefficient, reactants first;
//   - reactants and products are the two halves of compounds.
type reactionDecomposer struct {
	separatorPos int
	initCoefs    []float64
	compounds    []string
	reactants    []string
	products     []string
}

// newReactionDecomposer splits reaction on the "=" separator and then on "+",
// and separates the leading coefficient from each formula. It fails if the
// string is too short to hold a reaction, has no separator between reactants
// and products, or contains an empty compound (typically two adjacent "+").
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

// splitCoefFromFormula splits the leading coefficient (integer or float) off
// a compound of the reaction string and returns the pair (coefficient,
// formula). A compound that does not start with a digit keeps the coefficient 1.0.
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
