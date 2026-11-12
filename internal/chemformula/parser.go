package chemformula

import (
	"fmt"
	"math"
	"math/big"
	"slices"
	"strconv"
	"strings"
)

// Atom is a single entry of a parsed chemical formula: Label is the element
// symbol (e.g. "Fe") and Amount is the number of atoms of that element in the
// formula. Amount may be fractional, since the parser accepts non-integer
// atom counts such as 0.99 in solid-solution formulas.
type Atom struct {
	Label  string
	Amount float64
}

// String formats the atom as a dict-like entry, e.g. 'Fe': 2.
func (a Atom) String() string {
	return fmt.Sprintf("'%s': %v", a.Label, a.Amount)
}

// ExactAtom is a single entry of an exactly parsed chemical formula. Amount is
// an exact rational number of atoms, so nested fractional groups such as
// (Ga0.9Sc0.1)3 do not acquire binary floating-point artifacts.
type ExactAtom struct {
	Label  string
	Amount *big.Rat
}

// String formats the exact atom as a dict-like entry, e.g. 'Fe': 2.
func (a ExactAtom) String() string {
	return fmt.Sprintf("'%s': %s", a.Label, a.Amount.RatString())
}

// formulaCount carries the float and exact values of an atom count in
// parallel. The float value preserves the existing [Atom] API and behaviour;
// the rational value is used by the exact reaction balancer.
type formulaCount struct {
	float float64
	exact *big.Rat
}

// exactNumber parses a formula multiplier as an exact rational number.
func exactNumber(s string) (*big.Rat, error) {
	if strings.HasPrefix(s, ".") {
		s = "0" + s
	}
	r, ok := new(big.Rat).SetString(s)
	if !ok {
		return nil, fmt.Errorf("invalid number %q", s)
	}
	return r, nil
}

func formulaCountFrom(f float64, exact *big.Rat) formulaCount {
	if exact == nil {
		exact = new(big.Rat)
	}
	return formulaCount{float: f, exact: new(big.Rat).Set(exact)}
}

// add returns c + other in both representations.
func (c formulaCount) add(other formulaCount) formulaCount {
	exact := new(big.Rat)
	if c.exact != nil {
		exact.Add(exact, c.exact)
	}
	if other.exact != nil {
		exact.Add(exact, other.exact)
	}
	return formulaCount{float: c.float + other.float, exact: exact}
}

// mul returns c * other in both representations.
func (c formulaCount) mul(other formulaCount) formulaCount {
	exact := new(big.Rat)
	if c.exact != nil && other.exact != nil {
		exact.Mul(c.exact, other.exact)
	}
	return formulaCount{float: c.float * other.float, exact: exact}
}

// TokenType identifies the kind of lexical token produced by a [Lexer].
type TokenType int

const (
	// TokenElement is an element symbol such as "Fe" or "C".
	TokenElement TokenType = iota
	// TokenNumber is a numeric literal such as "2" or "0.5", used as an
	// atom amount or as the multiplier of a group or adduct.
	TokenNumber
	// TokenOpenParen is an opening bracket; after sanitization every
	// bracket style is reported as "(".
	TokenOpenParen
	// TokenCloseParen is a closing bracket; after sanitization every
	// bracket style is reported as ")".
	TokenCloseParen
	// TokenAdduct is the adduct separator: "*", "·" and "•" are all
	// normalized to "*", as in CuSO4*5H2O (water of crystallization
	// notation).
	TokenAdduct
	// TokenEOF marks the end of the input.
	TokenEOF
	// TokenInvalid is any character that cannot start a token, such as a
	// lowercase letter where an element symbol is expected.
	TokenInvalid
)

// Token is a lexical unit of a formula string. Position is the zero-based
// rune offset of the token in the normalized input; parser error messages
// report it as a 1-based position.
type Token struct {
	Type     TokenType
	Value    string
	Position int // Zero-based rune offset in the normalized input.
}

// Lexer splits a normalized formula string into [Token]s.
type Lexer struct {
	input []rune
	pos   int
}

// NewLexer returns a Lexer positioned at the start of input. Input must
// already be normalized by [sanitize].
func NewLexer(input []rune) *Lexer {
	return &Lexer{input: input, pos: 0}
}

// NextToken returns the next token of the input, or a [TokenEOF] token once
// the input is exhausted. Characters that cannot start a valid token are
// returned one by one as [TokenInvalid].
func (l *Lexer) NextToken() Token {
	if l.pos >= len(l.input) {
		return Token{Type: TokenEOF, Position: l.pos}
	}

	start := l.pos
	ch := l.input[l.pos]

	switch ch {
	case '(':
		l.pos++
		return Token{Type: TokenOpenParen, Value: "(", Position: start}
	case ')':
		l.pos++
		return Token{Type: TokenCloseParen, Value: ")", Position: start}
	case '*':
		l.pos++
		return Token{Type: TokenAdduct, Value: "*", Position: start}
	}
	if (ch >= '0' && ch <= '9') || ch == '.' {
		return l.readNumber()
	}

	if ch >= 'A' && ch <= 'Z' {
		return l.readElement()
	}

	l.pos++
	return Token{Type: TokenInvalid, Value: string(ch), Position: start}
}

// readElement reads an element symbol: an uppercase letter followed by any
// run of lowercase letters ("Uuo", "Fe").
func (l *Lexer) readElement() Token {
	start := l.pos
	l.pos++

	for l.pos < len(l.input) && l.input[l.pos] >= 'a' && l.input[l.pos] <= 'z' {
		l.pos++
	}

	return Token{Type: TokenElement, Value: string(l.input[start:l.pos]), Position: start}
}

// readNumber reads a numeric literal: a run of digits and decimal points.
// Whether the result is a valid number (at most one point, no trailing
// point) is decided later by [Parser.parseMultiplier].
func (l *Lexer) readNumber() Token {
	start := l.pos

	for l.pos < len(l.input) {
		ch := l.input[l.pos]
		if (ch >= '0' && ch <= '9') || ch == '.' {
			l.pos++
		} else {
			break
		}
	}

	return Token{Type: TokenNumber, Value: string(l.input[start:l.pos]), Position: start}
}

// Parser turns a normalized chemical formula string into a slice of [Atom]s.
// Amounts of atoms that appear several times are summed, parenthesized groups are
// folded into their surroundings with their multiplier, a single adduct ("*")
// joins two parts of the formula, and the resulting atoms keep the order of
// their first appearance in the string.
// It rejects empty formulas, unknown element symbols, stray characters,
// unbalanced parentheses, empty groups, malformed numbers and more than
// one adduct symbol. Errors reported while scanning carry a 1-based position
// within the normalized formula.
type Parser struct {
	lexer        *Lexer
	current      Token
	elementOrder []string
	seen         map[string]bool
}

// NewParser returns a Parser for formula, which must already be normalized
// by [sanitize]. The first token is read eagerly.
func NewParser(formula string) *Parser {
	lexer := NewLexer([]rune(formula))
	return &Parser{
		lexer:   lexer,
		current: lexer.NextToken(),
	}
}

// advance reads the next token from the lexer.
func (p *Parser) advance() {
	p.current = p.lexer.NextToken()
}

// parse parses the whole formula and returns its atoms ordered by first
// appearance, or an error describing the first problem it encountered. It
// backs the validation performed by [NewChemicalFormula].
func (p *Parser) parse() ([]Atom, error) {
	atoms, _, err := p.parseWithExact()
	return atoms, err
}

// parseWithExact parses the whole formula and returns both the float atoms used
// by the existing [Atom] API and the exact rational atoms used by the exact
// reaction balancer. The two representations are built in one pass so they
// cannot diverge.
func (p *Parser) parseWithExact() ([]Atom, []ExactAtom, error) {
	if p.current.Type == TokenEOF {
		return nil, nil, fmt.Errorf("empty formula string")
	}
	p.seen = make(map[string]bool)
	p.elementOrder = nil

	atomCounts, err := p.parseSequence(false)
	if err != nil {
		return nil, nil, err
	}

	if p.current.Type == TokenAdduct {
		adduct := p.current
		p.advance()
		multiplier, err := p.parseMultiplier()
		if err != nil {
			return nil, nil, err
		}
		subCounts, err := p.parseSequence(false)
		if err != nil {
			return nil, nil, err
		}
		if err := p.mergeCounts(atomCounts, subCounts, multiplier, adduct); err != nil {
			return nil, nil, err
		}
	}

	if p.current.Type == TokenAdduct {
		return nil, nil, p.errorAt(p.current, "there are more than 1 adduct symbol *•·")
	}
	if p.current.Type != TokenEOF {
		return nil, nil, p.unexpectedToken()
	}

	result := make([]Atom, 0, len(p.elementOrder))
	exactResult := make([]ExactAtom, 0, len(p.elementOrder))
	for _, label := range p.elementOrder {
		count := atomCounts[label]
		result = append(result, Atom{Label: label, Amount: count.float})
		exact := count.exact
		if exact == nil {
			exact = new(big.Rat)
		}
		exactResult = append(exactResult, ExactAtom{
			Label:  label,
			Amount: new(big.Rat).Set(exact),
		})
	}
	return result, exactResult, nil
}

// parseSequence parses atoms, parenthesized groups and their multipliers
// until a closing parenthesis (when inGroup is true), an adduct symbol or the
// end of input is reached, and returns the atom counts it accumulated. Each
// element symbol is recorded in p.elementOrder on first sight so that the
// final [Atom] slice preserves the order of the formula.
func (p *Parser) parseSequence(inGroup bool) (map[string]formulaCount, error) {
	atomCounts := make(map[string]formulaCount)
	terms := 0
	for {
		token := p.current
		switch token.Type {
		case TokenElement:
			if !isValidElement(token.Value) {
				return nil, p.errorAt(token, "there are invalid atom(s) %s", token.Value)
			}
			p.advance()
			count, err := p.parseMultiplier()
			if err != nil {
				return nil, err
			}
			if err := p.addCount(atomCounts, token.Value, count, token); err != nil {
				return nil, err
			}
			if !p.seen[token.Value] {
				p.elementOrder = append(p.elementOrder, token.Value)
				p.seen[token.Value] = true
			}

		case TokenOpenParen:
			p.advance()
			subCounts, err := p.parseSequence(true)
			if err != nil {
				return nil, err
			}
			// parseSequence(true) succeeds only at a closing parenthesis.
			p.advance()
			multiplier, err := p.parseMultiplier()
			if err != nil {
				return nil, err
			}
			if err := p.mergeCounts(atomCounts, subCounts, multiplier, token); err != nil {
				return nil, err
			}

		case TokenCloseParen:
			if !inGroup {
				return nil, p.errorAt(token, "parentheses [{()}] are not balanced: unexpected closing parenthesis")
			}
			if terms == 0 {
				return nil, p.errorAt(token, "empty parentheses group ()")
			}
			return atomCounts, nil

		case TokenEOF, TokenAdduct:
			if inGroup {
				if token.Type == TokenAdduct {
					return nil, p.errorAt(token, "adduct symbols are only allowed at the top level")
				}
				return nil, p.errorAt(token, "parentheses [{()}] are not balanced: missing closing parenthesis")
			}
			if terms == 0 {
				return nil, p.errorAt(token, "expected a nonempty formula sequence")
			}
			return atomCounts, nil

		default:
			return nil, p.unexpectedToken()
		}
		terms++
	}
}

// parseMultiplier consumes the numeric token that follows an atom, a group or
// an adduct symbol and returns its float and exact values, defaulting to 1
// when no number is present. Numbers with several decimal points, a trailing
// point or a non-finite float value are rejected.
func (p *Parser) parseMultiplier() (formulaCount, error) {
	if p.current.Type != TokenNumber {
		return formulaCountFrom(1, big.NewRat(1, 1)), nil
	}
	token := p.current
	if strings.Count(token.Value, ".") > 1 || strings.HasSuffix(token.Value, ".") {
		return formulaCount{}, p.errorAt(token, "invalid number %q", token.Value)
	}
	count, err := strconv.ParseFloat(token.Value, 64)
	if err != nil {
		return formulaCount{}, fmt.Errorf("%v: %w", p.errorAt(token, "invalid number %q", token.Value), err)
	}
	if math.IsNaN(count) || math.IsInf(count, 0) {
		return formulaCount{}, p.errorAt(token, "non-finite number %q", token.Value)
	}
	exact, err := exactNumber(token.Value)
	if err != nil {
		return formulaCount{}, fmt.Errorf("%v: %w", p.errorAt(token, "invalid number %q", token.Value), err)
	}
	p.advance()
	return formulaCountFrom(count, exact), nil
}

// addCount adds amount to the count of the element label, rejecting totals
// that would become NaN or infinite in the float representation.
func (p *Parser) addCount(counts map[string]formulaCount, label string, amount formulaCount, token Token) error {
	count := counts[label].add(amount)
	if math.IsNaN(count.float) || math.IsInf(count.float, 0) {
		return p.errorAt(token, "non-finite atom amount for %q", label)
	}
	counts[label] = count
	return nil
}

// mergeCounts folds the atom counts of a parsed group or adduct part into
// counts, scaling them by multiplier; token is used for error reporting.
func (p *Parser) mergeCounts(counts, subCounts map[string]formulaCount, multiplier formulaCount, token Token) error {
	for label, count := range subCounts {
		scaled := count.mul(multiplier)
		if math.IsNaN(scaled.float) || math.IsInf(scaled.float, 0) {
			return p.errorAt(token, "non-finite atom amount for %q", label)
		}
		if err := p.addCount(counts, label, scaled, token); err != nil {
			return err
		}
	}
	return nil
}

// unexpectedToken builds an error for a token that cannot appear at the
// current position, tailoring the message: a number in a place where an
// atom belongs ("No letters A-Z or a-z" if the formula contains no letters
// at all) or an invalid character (a lowercase letter is reported as an
// invalid atom).
func (p *Parser) unexpectedToken() error {
	token := p.current
	if token.Type == TokenNumber {
		if !slices.ContainsFunc(p.lexer.input, isLetter) {
			return p.errorAt(token, "no letters A-Z or a-z")
		}
		return p.errorAt(token, "unexpected number %q", token.Value)
	}
	if token.Type == TokenInvalid {
		if token.Value[0] >= 'a' && token.Value[0] <= 'z' {
			return p.errorAt(token, "there are invalid atom(s) %s", token.Value)
		}
		return p.errorAt(token, "there are invalid character(s) %s", token.Value)
	}
	return p.errorAt(token, "unexpected token %q", token.Value)
}

// errorAt wraps message with the 1-based position of token and the text of
// the normalized formula, e.g.:
// There are invalid atom(s) Xx at position 3 in formula "H2XxO".
func (p *Parser) errorAt(token Token, format string, args ...any) error {
	return fmt.Errorf("%s at position %d in formula %q", fmt.Sprintf(format, args...), token.Position+1, string(p.lexer.input))
}
