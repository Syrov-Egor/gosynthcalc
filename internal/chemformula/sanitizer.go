package chemformula

type formulaSanitizer struct{}

func (s formulaSanitizer) sanitize(formula string) string {
	res := make([]rune, 0, len(formula))
	for _, char := range formula {
		switch char {
		case ' ':
			continue
		case '[', '{':
			res = append(res, '(')
		case ']', '}':
			res = append(res, ')')
		case '·', '•':
			res = append(res, '*')
		default:
			res = append(res, char)
		}
	}
	return string(res)
}
