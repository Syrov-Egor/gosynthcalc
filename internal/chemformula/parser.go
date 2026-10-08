package chemformula

import (
	"fmt"
	"math"
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

// String formats the atom as a Python-dict-like entry, e.g. 'Fe': 2.
func (a Atom) String() string {
	return fmt.Sprintf("'%s': %v", a.Label, a.Amount)
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
	if p.current.Type == TokenEOF {
		return nil, fmt.Errorf("Empty formula string")
	}
	p.seen = make(map[string]bool)
	p.elementOrder = nil

	atomCounts, err := p.parseSequence(false)
	if err != nil {
		return nil, err
	}

	if p.current.Type == TokenAdduct {
		adduct := p.current
		p.advance()
		multiplier, err := p.parseMultiplier()
		if err != nil {
			return nil, err
		}
		subCounts, err := p.parseSequence(false)
		if err != nil {
			return nil, err
		}
		if err := p.mergeCounts(atomCounts, subCounts, multiplier, adduct); err != nil {
			return nil, err
		}
	}

	if p.current.Type == TokenAdduct {
		return nil, p.errorAt(p.current, "There are more than 1 adduct symbol *•·")
	}
	if p.current.Type != TokenEOF {
		return nil, p.unexpectedToken()
	}

	result := make([]Atom, 0, len(p.elementOrder))
	for _, label := range p.elementOrder {
		result = append(result, Atom{Label: label, Amount: atomCounts[label]})
	}
	return result, nil
}

// parseSequence parses atoms, parenthesized groups and their multipliers
// until a closing parenthesis (when inGroup is true), an adduct symbol or the
// end of input is reached, and returns the atom counts it accumulated. Each
// element symbol is recorded in p.elementOrder on first sight so that the
// final [Atom] slice preserves the order of the formula.
func (p *Parser) parseSequence(inGroup bool) (map[string]float64, error) {
	atomCounts := make(map[string]float64)
	terms := 0
	for {
		token := p.current
		switch token.Type {
		case TokenElement:
			if !isValidElement(token.Value) {
				return nil, p.errorAt(token, "There are invalid atom(s) %s", token.Value)
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
				return nil, p.errorAt(token, "Parentheses [{()}] are not balanced: unexpected closing parenthesis")
			}
			if terms == 0 {
				return nil, p.errorAt(token, "Empty parentheses group ()")
			}
			return atomCounts, nil

		case TokenEOF, TokenAdduct:
			if inGroup {
				if token.Type == TokenAdduct {
					return nil, p.errorAt(token, "Adduct symbols are only allowed at the top level")
				}
				return nil, p.errorAt(token, "Parentheses [{()}] are not balanced: missing closing parenthesis")
			}
			if terms == 0 {
				return nil, p.errorAt(token, "Expected a nonempty formula sequence")
			}
			return atomCounts, nil

		default:
			return nil, p.unexpectedToken()
		}
		terms++
	}
}

// parseMultiplier consumes the numeric token that follows an atom, a group or
// an adduct symbol and returns its value, defaulting to 1 when no number is
// present. Numbers with several decimal points, a trailing point or a
// non-finite value are rejected.
func (p *Parser) parseMultiplier() (float64, error) {
	if p.current.Type != TokenNumber {
		return 1, nil
	}
	token := p.current
	if strings.Count(token.Value, ".") > 1 || strings.HasSuffix(token.Value, ".") {
		return 0, p.errorAt(token, "invalid number %q", token.Value)
	}
	count, err := strconv.ParseFloat(token.Value, 64)
	if err != nil {
		return 0, fmt.Errorf("%v: %w", p.errorAt(token, "invalid number %q", token.Value), err)
	}
	if math.IsNaN(count) || math.IsInf(count, 0) {
		return 0, p.errorAt(token, "non-finite number %q", token.Value)
	}
	p.advance()
	return count, nil
}

// addCount adds amount to the count of the element label, rejecting totals
// that would become NaN or infinite.
func (p *Parser) addCount(counts map[string]float64, label string, amount float64, token Token) error {
	count := counts[label] + amount
	if math.IsNaN(count) || math.IsInf(count, 0) {
		return p.errorAt(token, "non-finite atom amount for %q", label)
	}
	counts[label] = count
	return nil
}

// mergeCounts folds the atom counts of a parsed group or adduct part into
// counts, scaling them by multiplier; token is used for error reporting.
func (p *Parser) mergeCounts(counts, subCounts map[string]float64, multiplier float64, token Token) error {
	for label, count := range subCounts {
		if err := p.addCount(counts, label, count*multiplier, token); err != nil {
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
			return p.errorAt(token, "No letters A-Z or a-z")
		}
		return p.errorAt(token, "Unexpected number %q", token.Value)
	}
	if token.Type == TokenInvalid {
		if token.Value[0] >= 'a' && token.Value[0] <= 'z' {
			return p.errorAt(token, "There are invalid atom(s) %s", token.Value)
		}
		return p.errorAt(token, "There are invalid character(s) %s", token.Value)
	}
	return p.errorAt(token, "Unexpected token %q", token.Value)
}

// errorAt wraps message with the 1-based position of token and the text of
// the normalized formula, e.g.:
// There are invalid atom(s) Xx at position 3 in formula "H2XxO".
func (p *Parser) errorAt(token Token, format string, args ...any) error {
	return fmt.Errorf("%s at position %d in formula %q", fmt.Sprintf(format, args...), token.Position+1, string(p.lexer.input))
}
