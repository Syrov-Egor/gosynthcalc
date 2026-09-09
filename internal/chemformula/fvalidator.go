package chemformula

import (
	"fmt"
	"slices"
	"strings"
)

var (
	allowed     [256]bool
	validSingle [26]bool
	validDouble [26][26]bool
)

func init() {
	for c := 'a'; c <= 'z'; c++ {
		allowed[c] = true
	}
	for c := 'A'; c <= 'Z'; c++ {
		allowed[c] = true
	}
	for c := '0'; c <= '9'; c++ {
		allowed[c] = true
	}
	for _, c := range "().*" {
		allowed[c] = true
	}

	for k := range periodicTable {
		switch len(k) {
		case 1:
			validSingle[k[0]-'A'] = true
		case 2:
			validDouble[k[0]-'A'][k[1]-'a'] = true
		}
	}
}

func sanitize(formula string) string {
	var res strings.Builder
	for _, r := range formula {
		switch r {
		case '[', '{':
			res.WriteRune('(')
		case ']', '}':
			res.WriteRune(')')
		case '·', '•':
			res.WriteRune('*')
		case ',':
			res.WriteRune('.')
		case ' ':
		default:
			res.WriteRune(r)
		}
	}
	return res.String()
}

type formulaValidator struct {
	formula string
}

func (v formulaValidator) validate() error {

	if v.formula == "" {
		return fmt.Errorf("Empty formula string")
	}

	leftParenthesisCount, rightParenthesisCount, adductCount := 0, 0, 0
	letterPresent := false
	invalidCharacters := make([]rune, 0)
	var allLetters strings.Builder

	for _, r := range v.formula {
		if !letterPresent {
			letterPresent = isLetter(r)
		}
		if isLetter(r) {
			allLetters.WriteRune(r)
		}
		if !v.isAllowed(r) {
			invalidCharacters = append(invalidCharacters, r)
		}
		switch r {
		case '(':
			leftParenthesisCount++
		case ')':
			rightParenthesisCount++
		case '*':
			adductCount++
		}
	}

	if !letterPresent {
		return fmt.Errorf("No letters A-Z or a-z in the formula '%s'", v.formula)
	}
	if len(invalidCharacters) > 0 {
		return fmt.Errorf("There are invalid character(s) %s in the formula '%s'", string(invalidCharacters), v.formula)
	}
	invalidAtoms := invalidAtoms(v.formula)
	if len(invalidAtoms) > 0 {
		return fmt.Errorf("There are invalid atom(s) %s in the formula '%s'", invalidAtoms, v.formula)
	}
	if leftParenthesisCount != rightParenthesisCount {
		return fmt.Errorf("Parentheses [{()}] are not balanced in the formula '%s'", v.formula)
	}
	if adductCount > 1 {
		return fmt.Errorf("There are more than 1 adduct symbol *•· in the formula '%s'", v.formula)
	}
	return nil
}

func (v formulaValidator) isAllowed(r rune) bool {
	return r >= 0 && r < 256 && allowed[byte(r)]
}

func isLetter(r rune) bool {
	return (r >= 'a' && r <= 'z') || (r >= 'A' && r <= 'Z')
}

func isValidElement(tok string) bool {
	switch len(tok) {
	case 1:
		c := tok[0]
		return c >= 'A' && c <= 'Z' && validSingle[c-'A']
	case 2:
		a, b := tok[0], tok[1]
		return a >= 'A' && a <= 'Z' && b >= 'a' && b <= 'z' &&
			validDouble[a-'A'][b-'a']
	default:
		return false
	}
}

func invalidAtoms(text string) []string {
	var res []string
	for i := 0; i < len(text); {
		c := text[i]
		if c >= 'A' && c <= 'Z' {
			j := i + 1
			for j < len(text) && text[j] >= 'a' && text[j] <= 'z' {
				j++
			}
			tok := text[i:j]
			if !isValidElement(tok) {
				dup := slices.Contains(res, tok)
				if !dup {
					res = append(res, tok)
				}
			}
			i = j
		} else if c >= 'a' && c <= 'z' {
			tok := text[i : i+1]
			dup := slices.Contains(res, tok)
			if !dup {
				res = append(res, tok)
			}
			i++
		} else {
			i++
		}
	}
	return res
}
