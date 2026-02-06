package chemformula

import (
	"strconv"
	"unicode"
)

type TokenType int

const (
	TokenElement TokenType = iota
	TokenNumber
	TokenOpenParen
	TokenCloseParen
	TokenAdduct
	TokenEOF
)

type Token struct {
	Type  TokenType
	Value string
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
		return Token{Type: TokenEOF}
	}

	ch := l.input[l.pos]

	switch ch {
	case '(':
		l.pos++
		return Token{Type: TokenOpenParen, Value: "("}
	case ')':
		l.pos++
		return Token{Type: TokenCloseParen, Value: ")"}
	case '*':
		l.pos++
		return Token{Type: TokenAdduct, Value: "*"}
	}
	if unicode.IsDigit(ch) || ch == '.' {
		return l.readNumber()
	}

	if unicode.IsUpper(ch) {
		return l.readElement()
	}

	l.pos++
	return l.NextToken()
}

func (l *Lexer) readElement() Token {
	start := l.pos
	l.pos++

	for l.pos < len(l.input) && unicode.IsLower(rune(l.input[l.pos])) {
		l.pos++
	}

	return Token{Type: TokenElement, Value: string(l.input[start:l.pos])}
}

func (l *Lexer) readNumber() Token {
	start := l.pos

	for l.pos < len(l.input) {
		ch := rune(l.input[l.pos])
		if unicode.IsDigit(ch) || ch == '.' {
			l.pos++
		} else {
			break
		}
	}

	return Token{Type: TokenNumber, Value: string(l.input[start:l.pos])}
}

type Parser struct {
	lexer   *Lexer
	current Token
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

func (p *Parser) parse() []Atom {
	atomCounts := make(map[string]float64, len(p.lexer.input)/2)
	elementOrder := []string{}
	seen := make(map[string]bool)

	p.parseFormula(atomCounts, &elementOrder, &seen, 1.0)

	var result []Atom
	for _, label := range elementOrder {
		if count, exists := atomCounts[label]; exists {
			result = append(result, Atom{Label: label, Amount: count})
		}
	}

	return result
}

func (p *Parser) parseFormula(atomCounts map[string]float64,
	elementOrder *[]string,
	seen *map[string]bool,
	multiplier float64) {
	for p.current.Type != TokenEOF {
		switch p.current.Type {
		case TokenElement:
			element := p.current.Value
			p.advance()

			count := 1.0
			if p.current.Type == TokenNumber {
				count, _ = strconv.ParseFloat(p.current.Value, 64)
				p.advance()
			}

			atomCounts[element] += count * multiplier

			if !(*seen)[element] {
				*elementOrder = append(*elementOrder, element)
				(*seen)[element] = true
			}

		case TokenOpenParen:
			p.advance()
			subCounts := make(map[string]float64)
			p.parseGroup(subCounts, elementOrder, seen, 1.0)

			groupMultiplier := 1.0
			if p.current.Type == TokenNumber {
				groupMultiplier, _ = strconv.ParseFloat(p.current.Value, 64)
				p.advance()
			}

			for label, count := range subCounts {
				atomCounts[label] += count * groupMultiplier * multiplier
			}

		case TokenAdduct:
			p.advance()

			adductMultiplier := 1.0
			if p.current.Type == TokenNumber {
				adductMultiplier, _ = strconv.ParseFloat(p.current.Value, 64)
				p.advance()
			}

			p.parseFormula(atomCounts, elementOrder, seen, adductMultiplier)

		default:
			p.advance()
		}
	}
}

func (p *Parser) parseGroup(atomCounts map[string]float64, elementOrder *[]string, seen *map[string]bool, multiplier float64) {
	for p.current.Type != TokenEOF {
		switch p.current.Type {
		case TokenCloseParen:
			p.advance()
			return

		case TokenElement:
			element := p.current.Value
			p.advance()

			count := 1.0
			if p.current.Type == TokenNumber {
				count, _ = strconv.ParseFloat(p.current.Value, 64)
				p.advance()
			}

			atomCounts[element] += count * multiplier

			if !(*seen)[element] {
				*elementOrder = append(*elementOrder, element)
				(*seen)[element] = true
			}

		case TokenOpenParen:
			p.advance()
			subCounts := make(map[string]float64)
			p.parseGroup(subCounts, elementOrder, seen, 1.0)

			groupMultiplier := 1.0
			if p.current.Type == TokenNumber {
				groupMultiplier, _ = strconv.ParseFloat(p.current.Value, 64)
				p.advance()
			}

			for label, count := range subCounts {
				atomCounts[label] += count * groupMultiplier * multiplier
			}

		default:
			p.advance()
		}
	}
}
