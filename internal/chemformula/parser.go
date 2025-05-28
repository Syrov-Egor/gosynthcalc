package chemformula

import (
	"fmt"
	"regexp"
	"slices"
	"strconv"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

type regexes struct {
	atomRegex        *regexp.Regexp
	coefRegex        *regexp.Regexp
	atomAndCoefRegex *regexp.Regexp
	letterRegex      *regexp.Regexp
	noLetterRegex    *regexp.Regexp
	allowedSymbols   *regexp.Regexp
	openerBrackets   []rune
	closerBrackets   []rune
	adductSymbols    []rune
}

var formRegexes regexes = regexes{
	atomRegex:        regexp.MustCompile(`([A-Z][a-z]*)`),
	coefRegex:        regexp.MustCompile(`((\d+(\.\d+)?)*)`),
	atomAndCoefRegex: regexp.MustCompile(`([A-Z][a-z]*)((\d+(\.\d+)?)*)`),
	letterRegex:      regexp.MustCompile(`[a-z]`),
	noLetterRegex:    regexp.MustCompile(`[A-Za-z]`),
	allowedSymbols:   regexp.MustCompile(`[^A-Za-z0-9.({[)}\]*·•]`),
	openerBrackets:   []rune{'(', '[', '{'},
	closerBrackets:   []rune{')', ']', '}'},
	adductSymbols:    []rune{'*', '·', '•'},
}

type Atom struct {
	Label  string
	Amount float64
}

func (a Atom) String() string {
	return fmt.Sprintf("'%s': %v", a.Label, a.Amount)
}

type chemicalFormulaParser struct{}

func (p chemicalFormulaParser) parseToMap(formula string) (map[string]float64, int) {
	tokens := []rune{}
	mol := make(map[string]float64)
	i := 0

	for i < len(formula) {
		token := rune(formula[i])
		switch {

		case slices.Contains(formRegexes.adductSymbols, token):
			matches := formRegexes.coefRegex.FindStringSubmatch(formula[i+1:])
			weight := 1.0

			if len(matches) > 0 && matches[0] != "" {
				weight, _ = strconv.ParseFloat(matches[0], 32)
				i += len(matches[0])
			}

			submol, lenght := p.parseToMap("(" + formula[i+1:] + ")" + strconv.FormatFloat(weight, 'f', -1, 64))
			mol = p.fuse(mol, submol, 1.0)
			i += lenght + 1

		case slices.Contains(formRegexes.closerBrackets, token):
			matches := formRegexes.coefRegex.FindStringSubmatch(formula[i+1:])
			weight := 1.0

			if len(matches) > 0 && matches[0] != "" {
				weight, _ = strconv.ParseFloat(matches[0], 64)
				i += len(matches[0])
			}

			tokenStr := string(tokens)
			submol := p.toMap(formRegexes.atomAndCoefRegex.FindAllStringSubmatch(tokenStr, -1))
			return p.fuse(mol, submol, weight), i

		case slices.Contains(formRegexes.openerBrackets, token):
			submol, length := p.parseToMap(formula[i+1:])
			mol = p.fuse(mol, submol, 1.0)
			i += length + 1

		default:
			tokens = append(tokens, token)
		}
		i++
	}
	tokenStr := string(tokens)
	extractFromTokens := formRegexes.atomAndCoefRegex.FindAllStringSubmatch(tokenStr, -1)
	fusedMap := p.fuse(mol, p.toMap(extractFromTokens), 1.0)

	return fusedMap, i
}

func (p chemicalFormulaParser) fuse(mol1, mol2 map[string]float64, weight float64) map[string]float64 {
	fused := make(map[string]float64)
	for atom, count := range mol1 {
		fused[atom] += count * weight
	}
	for atom, count := range mol2 {
		fused[atom] += count * weight
	}
	return fused
}

func (p chemicalFormulaParser) toMap(matches [][]string) map[string]float64 {
	result := make(map[string]float64)

	for _, match := range matches {
		atom := match[1]
		nStr := match[2]
		var n float64 = 1.0
		if nStr != "" {
			var err error
			n, err = strconv.ParseFloat(nStr, 64)
			if err != nil {
				n = 1.0
			}
		}
		result[atom] += n
	}

	return result
}

func (p chemicalFormulaParser) order(formula string, parsed map[string]float64) []Atom {
	ret := make([]Atom, len(parsed))
	atomMatch := formRegexes.atomRegex.FindAllString(formula, -1)
	unique := utils.UniqueElems(atomMatch)
	for i, match := range unique {
		ret[i] = Atom{Label: match, Amount: parsed[match]}
	}
	return ret
}

func (p chemicalFormulaParser) parse(formula string) []Atom {
	parsed, _ := p.parseToMap(formula)
	res := p.order(formula, parsed)
	return res
}
