package chemformula

import (
	"fmt"
	"math"
	"slices"
	"strconv"
	"strings"
)

type Atom struct {
	Label  string
	Amount float64
}

func (a Atom) String() string {
	return fmt.Sprintf("'%s': %v", a.Label, a.Amount)
}

type TokenType int

const (
	TokenElement TokenType = iota
	TokenNumber
	TokenOpenParen
	TokenCloseParen
	TokenAdduct
	TokenEOF
	TokenInvalid
)

type Token struct {
	Type     TokenType
	Value    string
	Position int // Zero-based rune offset in the normalized input.
}

type Lexer struct {
	input []rune
	pos   int
}

func NewLexer(input []rune) *Lexer {
	return &Lexer{input: input, pos: 0}
}

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

func (l *Lexer) readElement() Token {
	start := l.pos
	l.pos++

	for l.pos < len(l.input) && l.input[l.pos] >= 'a' && l.input[l.pos] <= 'z' {
		l.pos++
	}

	return Token{Type: TokenElement, Value: string(l.input[start:l.pos]), Position: start}
}

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

type Parser struct {
	lexer        *Lexer
	current      Token
	elementOrder []string
	seen         map[string]bool
}

func NewParser(formula string) *Parser {
	lexer := NewLexer([]rune(formula))
	return &Parser{
		lexer:   lexer,
		current: lexer.NextToken(),
	}
}

func (p *Parser) advance() {
	p.current = p.lexer.NextToken()
}

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

func (p *Parser) addCount(counts map[string]float64, label string, amount float64, token Token) error {
	count := counts[label] + amount
	if math.IsNaN(count) || math.IsInf(count, 0) {
		return p.errorAt(token, "non-finite atom amount for %q", label)
	}
	counts[label] = count
	return nil
}

func (p *Parser) mergeCounts(counts, subCounts map[string]float64, multiplier float64, token Token) error {
	for label, count := range subCounts {
		if err := p.addCount(counts, label, count*multiplier, token); err != nil {
			return err
		}
	}
	return nil
}

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

func (p *Parser) errorAt(token Token, format string, args ...any) error {
	return fmt.Errorf("%s at position %d in formula %q", fmt.Sprintf(format, args...), token.Position+1, string(p.lexer.input))
}
