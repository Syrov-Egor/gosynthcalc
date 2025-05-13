package formula

import (
	"fmt"
	"regexp"
	"slices"
	"strconv"
)

var atomRegex *regexp.Regexp = regexp.MustCompile(`([A-Z][a-z]*)`)
var coefRegex *regexp.Regexp = regexp.MustCompile(`((\d+(\.\d+)?)*)`)
var atomAndCoefRegex *regexp.Regexp = regexp.MustCompile(`([A-Z][a-z]*)((\d+(\.\d+)?)*)`)
var openerBrackets []rune = []rune{'(', '[', '{'}
var closerBrackets []rune = []rune{')', ']', '}'}
var adductSymbols []rune = []rune{'*', '·', '•'}

func Parse(formula string) (map[string]float64, int) {
	tokens := []rune{}
	mol := make(map[string]float64)
	i := 0

	for i < len(formula) {
		token := rune(formula[i])
		switch {

		case slices.Contains(adductSymbols, token):
			matches := coefRegex.FindStringSubmatch(formula[i+1:])
			weight := 1.0

			if len(matches) > 0 && matches[0] != "" {
				weight, _ = strconv.ParseFloat(matches[0], 32)
				i += len(matches[0])
			}

			submol, lenght := Parse("(" + formula[i+1:] + ")" + strconv.FormatFloat(weight, 'f', -1, 64))
			mol = fuse(mol, submol, 1.0)
			fmt.Println("Fused on adduct:", mol)
			i += lenght + 1

		case slices.Contains(closerBrackets, token):
			matches := coefRegex.FindStringSubmatch(formula[i+1:])
			weight := 1.0

			if len(matches) > 0 && matches[0] != "" {
				weight, _ = strconv.ParseFloat(matches[0], 64)
				i += len(matches[0])
			}

			tokenStr := string(tokens)
			submol := toMap(atomAndCoefRegex.FindAllStringSubmatch(tokenStr, -1))
			fmt.Println("mol:", mol)
			fmt.Println("submol:", submol)
			fmt.Println("Weight:", weight)
			fmt.Println("Fused on closer:", fuse(mol, submol, weight))
			return fuse(mol, submol, weight), i

		case slices.Contains(openerBrackets, token):
			submol, length := Parse(formula[i+1:])
			mol = fuse(mol, submol, 1.0)
			fmt.Println("Fused on opener:", mol)
			i += length + 1

		default:
			tokens = append(tokens, token)
		}
		i++
	}
	tokenStr := string(tokens)
	extractFromTokens := atomAndCoefRegex.FindAllStringSubmatch(tokenStr, -1)
	fusedMap := fuse(mol, toMap(extractFromTokens), 1.0)
	fmt.Println("Final fuse:", fusedMap)
	return fusedMap, i
}

func fuse(mol1, mol2 map[string]float64, weight float64) map[string]float64 {
	fused := make(map[string]float64)
	for atom, count := range mol1 {
		fused[atom] += count * weight
	}
	for atom, count := range mol2 {
		fused[atom] += count * weight
	}
	return fused
}

func toMap(matches [][]string) map[string]float64 {
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
