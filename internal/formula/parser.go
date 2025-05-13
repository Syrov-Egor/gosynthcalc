package formula

import (
	"regexp"
	"slices"
	"strconv"
)

type ChemicalFormulaParser struct {
	atomRegex        *regexp.Regexp
	coefRegex        *regexp.Regexp
	atomAndCoefRegex *regexp.Regexp
	openerBrackets   []rune
	closerBrackets   []rune
	adductSymbols    []rune
}

func NewChemicalFormulaParser() *ChemicalFormulaParser {
	return &ChemicalFormulaParser{
		atomRegex:        regexp.MustCompile(`([A-Z][a-z]*)`),
		coefRegex:        regexp.MustCompile(`((\d+(\.\d+)?)*)`),
		atomAndCoefRegex: regexp.MustCompile(`([A-Z][a-z]*)((\d+(\.\d+)?)*)`),
		openerBrackets:   []rune{'(', '[', '{'},
		closerBrackets:   []rune{')', ']', '}'},
		adductSymbols:    []rune{'*', '·', '•'},
	}
}

func (p *ChemicalFormulaParser) parse(formula string) (map[string]float64, int) {
	tokens := []rune{}
	mol := make(map[string]float64)
	i := 0

	for i < len(formula) {
		token := rune(formula[i])
		switch {

		case slices.Contains(p.adductSymbols, token):
			matches := p.coefRegex.FindStringSubmatch(formula[i+1:])
			weight := 1.0

			if len(matches) > 0 && matches[0] != "" {
				weight, _ = strconv.ParseFloat(matches[0], 32)
				i += len(matches[0])
			}

			submol, lenght := p.parse("(" + formula[i+1:] + ")" + strconv.FormatFloat(weight, 'f', -1, 64))
			mol = p.fuse(mol, submol, 1.0)
			i += lenght + 1

		case slices.Contains(p.closerBrackets, token):
			matches := p.coefRegex.FindStringSubmatch(formula[i+1:])
			weight := 1.0

			if len(matches) > 0 && matches[0] != "" {
				weight, _ = strconv.ParseFloat(matches[0], 64)
				i += len(matches[0])
			}

			tokenStr := string(tokens)
			submol := p.toMap(p.atomAndCoefRegex.FindAllStringSubmatch(tokenStr, -1))
			return p.fuse(mol, submol, weight), i

		case slices.Contains(p.openerBrackets, token):
			submol, length := p.parse(formula[i+1:])
			mol = p.fuse(mol, submol, 1.0)
			i += length + 1

		default:
			tokens = append(tokens, token)
		}
		i++
	}
	tokenStr := string(tokens)
	extractFromTokens := p.atomAndCoefRegex.FindAllStringSubmatch(tokenStr, -1)
	fusedMap := p.fuse(mol, p.toMap(extractFromTokens), 1.0)
	return fusedMap, i
}

func (p *ChemicalFormulaParser) fuse(mol1, mol2 map[string]float64, weight float64) map[string]float64 {
	fused := make(map[string]float64)
	for atom, count := range mol1 {
		fused[atom] += count * weight
	}
	for atom, count := range mol2 {
		fused[atom] += count * weight
	}
	return fused
}

func (p *ChemicalFormulaParser) toMap(matches [][]string) map[string]float64 {
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
