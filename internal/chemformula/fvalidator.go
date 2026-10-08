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
	res.Grow(len(formula))
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
	_, err := NewParser(v.formula).parse()
	return err
}

func validateBrackets(formula string) error {
	var stack []rune
	position := 0
	for _, r := range formula {
		position++
		switch r {
		case '(', '[', '{':
			stack = append(stack, r)
		case ')', ']', '}':
			if len(stack) == 0 {
				return fmt.Errorf("Parentheses [{()}] are not balanced: unexpected %q at position %d in formula %q", r, position, formula)
			}
			open := stack[len(stack)-1]
			if (r == ')' && open != '(') || (r == ']' && open != '[') || (r == '}' && open != '{') {
				return fmt.Errorf("Parentheses [{()}] are not balanced: mismatched %q and %q at position %d in formula %q", open, r, position, formula)
			}
			stack = stack[:len(stack)-1]
		}
	}
	if len(stack) > 0 {
		return fmt.Errorf("Parentheses [{()}] are not balanced: missing closing bracket in formula %q", formula)
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
	var result []string
	for i := 0; i < len(text); {
		char := text[i]
		if char >= 'A' && char <= 'Z' {
			j := i + 1
			for j < len(text) && text[j] >= 'a' && text[j] <= 'z' {
				j++
			}
			token := text[i:j]
			if !isValidElement(token) {
				duplicate := slices.Contains(result, token)
				if !duplicate {
					result = append(result, token)
				}
			}
			i = j
		} else if char >= 'a' && char <= 'z' {
			token := text[i : i+1]
			duplicate := slices.Contains(result, token)
			if !duplicate {
				result = append(result, token)
			}
			i++
		} else {
			i++
		}
	}
	return result
}
