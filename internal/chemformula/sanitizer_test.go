package chemformula

import (
	"testing"
)

func Test_formulaSanitizer_sanitize(t *testing.T) {
	tests := []struct {
		name    string
		formula string
		want    string
	}{
		{
			name:    "sanitize string",
			formula: "{K2}  2Mg2[(SO4)3Ho]2",
			want:    "(K2)2Mg2((SO4)3Ho)2",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			var s formulaSanitizer
			got := s.sanitize(tt.formula)
			if tt.want != got {
				t.Errorf("sanitize() = %v, want %v", got, tt.want)
			}
		})
	}
}
