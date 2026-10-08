package chemformula

import (
	"fmt"
	"slices"
	"strings"
)

// allowed holds the ASCII characters permitted in a formula: letters, digits
// and the punctuation ". ( ) *". validSingle and validDouble hold the
// one- and two-letter symbols of the periodic table. All three tables are
// built once in init: the element tables directly from periodicTable, so
// they can never drift out of sync with it.
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

// sanitize normalizes a formula string before parsing: whitespace is dropped,
// square and curly brackets become parentheses, the adduct symbols "·" and
// "•" become "*" and commas become decimal points. Everything else is kept
// as-is, so "K4[Fe(CN)6]·3H2O" and "CuSO4•5H2O" become ordinary
// parenthesis/adduct notation for the parser.
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

// formulaValidator checks a raw formula string. Almost every
// check is carried out by the [Parser] itself, so validate only
// has to run a full parse and report its first error.
type formulaValidator struct {
	formula string
}

// validate parses the formula from scratch and returns the first parsing
// error, if any.
func (v formulaValidator) validate() error {
	_, err := NewParser(v.formula).parse()
	return err
}

// validateBrackets reports whether every "(", "[" and "{" in formula is
// closed by the matching bracket type in the correct nesting order, e.g.
// "([)]" is rejected. It runs before [sanitize], so all three bracket styles
// are checked, and the reported positions are 1-based rune offsets.
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
				return fmt.Errorf("parentheses [{()}] are not balanced: unexpected %q at position %d in formula %q", r, position, formula)
			}
			open := stack[len(stack)-1]
			if (r == ')' && open != '(') || (r == ']' && open != '[') || (r == '}' && open != '{') {
				return fmt.Errorf("parentheses [{()}] are not balanced: mismatched %q and %q at position %d in formula %q", open, r, position, formula)
			}
			stack = stack[:len(stack)-1]
		}
	}
	if len(stack) > 0 {
		return fmt.Errorf("parentheses [{()}] are not balanced: missing closing bracket in formula %q", formula)
	}
	return nil
}

// isAllowed reports whether r is one of the ASCII characters a formula may
// contain: a letter, a digit or one of ". ( ) *".
func (v formulaValidator) isAllowed(r rune) bool {
	return r >= 0 && r < 256 && allowed[byte(r)]
}

// isLetter reports whether r is an ASCII letter.
func isLetter(r rune) bool {
	return (r >= 'a' && r <= 'z') || (r >= 'A' && r <= 'Z')
}

// isValidElement reports whether tok is a symbol known to the periodic
// table. Only one-letter and two-letter symbols in canonical case ("C",
// "Fe") are accepted, so "FE" or "x" are not.
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

// invalidAtoms returns the distinct element-like tokens of text that are not
// in the periodic table, each reported once, together with any stray
// lowercase letters. The string is scanned left to right for
// runs that start with an uppercase letter and continue with lowercase ones
// ("Zz" in "CO2Zz"), and letters that cannot form such a run are reported
// individually.
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
