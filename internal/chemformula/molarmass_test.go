package chemformula

import (
	"slices"
	"testing"

	"github.com/Syrov-Egor/gosynthcalc/internal/utils"
)

func TestMolarMass_molarMass(t *testing.T) {
	tests := []struct {
		name     string
		parsed   []Atom
		expected float64
	}{
		{
			name: "formula 1",
			parsed: []Atom{
				{Label: "H", Amount: 2},
				{Label: "O", Amount: 1}},
			expected: 18.015,
		},
		{
			name: "formula 2",
			parsed: []Atom{
				{Label: "N", Amount: 2},
				{Label: "H", Amount: 10},
				{Label: "S", Amount: 1},
				{Label: "O", Amount: 5}},
			expected: 150.149,
		},
		{
			name: "formula 3",
			parsed: []Atom{
				{Label: "K", Amount: 1.2},
				{Label: "Na", Amount: 0.8},
				{Label: "S", Amount: 1},
				{Label: "O", Amount: 4}},
			expected: 161.365415424,
		},
		{
			name: "formula 4",
			parsed: []Atom{
				{Label: "K", Amount: 4},
				{Label: "Mg", Amount: 2},
				{Label: "S", Amount: 6},
				{Label: "O", Amount: 24},
				{Label: "Ho", Amount: 2}},
			expected: 1111.198658,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			m := molarMass{tt.parsed}
			result := utils.RoundFloat(m.molarMass(), 10)
			if result != tt.expected {
				t.Errorf("molarMass() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestMolarMass_massPercent(t *testing.T) {
	tests := []struct {
		name     string
		parsed   []Atom
		expected []Atom
	}{
		{
			name: "formula 1",
			parsed: []Atom{
				{Label: "H", Amount: 2},
				{Label: "O", Amount: 1}},
			expected: []Atom{
				{Label: "H", Amount: 11.19067443796836},
				{Label: "O", Amount: 88.80932556203163}},
		},
		{
			name: "formula 2",
			parsed: []Atom{
				{Label: "N", Amount: 2},
				{Label: "H", Amount: 10},
				{Label: "S", Amount: 1},
				{Label: "O", Amount: 5}},
			expected: []Atom{
				{Label: "N", Amount: 18.657466916196576},
				{Label: "H", Amount: 6.713331424118708},
				{Label: "S", Amount: 21.35212355726645},
				{Label: "O", Amount: 53.277078102418265}},
		},
		{
			name: "formula 3",
			parsed: []Atom{
				{Label: "K", Amount: 1.2},
				{Label: "Na", Amount: 0.8},
				{Label: "S", Amount: 1},
				{Label: "O", Amount: 4}},
			expected: []Atom{
				{Label: "K", Amount: 29.075375213902188},
				{Label: "Na", Amount: 11.397619109196413},
				{Label: "S", Amount: 19.867949966701286},
				{Label: "O", Amount: 39.65905571020011}},
		},
		{
			name: "formula 4",
			parsed: []Atom{
				{Label: "K", Amount: 4},
				{Label: "Mg", Amount: 2},
				{Label: "S", Amount: 6},
				{Label: "O", Amount: 24},
				{Label: "Ho", Amount: 2}},
			expected: []Atom{
				{Label: "K", Amount: 14.074171065098604},
				{Label: "Mg", Amount: 4.374555319162381},
				{Label: "S", Amount: 17.31103602538728},
				{Label: "O", Amount: 34.5551173262846},
				{Label: "Ho", Amount: 29.685120264067127}},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			m := molarMass{tt.parsed}
			result := m.massPercent()
			if !slices.Equal(result, tt.expected) {
				t.Errorf("molarMass() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestMolarMass_atomicPercent(t *testing.T) {
	tests := []struct {
		name     string
		parsed   []Atom
		expected []Atom
	}{
		{
			name: "formula 1",
			parsed: []Atom{
				{Label: "H", Amount: 2},
				{Label: "O", Amount: 1}},
			expected: []Atom{
				{Label: "H", Amount: 66.66666666666666},
				{Label: "O", Amount: 33.33333333333333}},
		},
		{
			name: "formula 2",
			parsed: []Atom{
				{Label: "N", Amount: 2},
				{Label: "H", Amount: 10},
				{Label: "S", Amount: 1},
				{Label: "O", Amount: 5}},
			expected: []Atom{
				{Label: "N", Amount: 11.11111111111111},
				{Label: "H", Amount: 55.55555555555556},
				{Label: "S", Amount: 5.555555555555555},
				{Label: "O", Amount: 27.77777777777778}},
		},
		{
			name: "formula 3",
			parsed: []Atom{
				{Label: "K", Amount: 1.2},
				{Label: "Na", Amount: 0.8},
				{Label: "S", Amount: 1},
				{Label: "O", Amount: 4}},
			expected: []Atom{
				{Label: "K", Amount: 17.142857142857142},
				{Label: "Na", Amount: 11.428571428571429},
				{Label: "S", Amount: 14.285714285714285},
				{Label: "O", Amount: 57.14285714285714}},
		},
		{
			name: "formula 4",
			parsed: []Atom{
				{Label: "K", Amount: 4},
				{Label: "Mg", Amount: 2},
				{Label: "S", Amount: 6},
				{Label: "O", Amount: 24},
				{Label: "Ho", Amount: 2}},
			expected: []Atom{
				{Label: "K", Amount: 10.526315789473683},
				{Label: "Mg", Amount: 5.263157894736842},
				{Label: "S", Amount: 15.789473684210526},
				{Label: "O", Amount: 63.1578947368421},
				{Label: "Ho", Amount: 5.263157894736842}},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			m := molarMass{tt.parsed}
			result := m.atomicPercent()
			if !slices.Equal(result, tt.expected) {
				t.Errorf("molarMass() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestMolarMass_oxidePercent(t *testing.T) {
	tests := []struct {
		name     string
		parsed   []Atom
		expected []Atom
	}{
		{
			name: "formula 1",
			parsed: []Atom{
				{Label: "H", Amount: 2},
				{Label: "O", Amount: 1}},
			expected: []Atom{
				{Label: "H2O", Amount: 100}},
		},
		{
			name: "formula 2",
			parsed: []Atom{
				{Label: "N", Amount: 2},
				{Label: "H", Amount: 10},
				{Label: "S", Amount: 1},
				{Label: "O", Amount: 5}},
			expected: []Atom{
				{Label: "NO2", Amount: 35.09929732740271},
				{Label: "H2O", Amount: 34.36114777487011},
				{Label: "SO3", Amount: 30.53955489772719}},
		},
		{
			name: "formula 3",
			parsed: []Atom{
				{Label: "K", Amount: 1.2},
				{Label: "Na", Amount: 0.8},
				{Label: "S", Amount: 1},
				{Label: "O", Amount: 4}},
			expected: []Atom{
				{Label: "K2O", Amount: 35.02423357043221},
				{Label: "Na2O", Amount: 15.363524680216425},
				{Label: "SO3", Amount: 49.61224174935137}},
		},
		{
			name: "formula 4",
			parsed: []Atom{
				{Label: "K", Amount: 4},
				{Label: "Mg", Amount: 2},
				{Label: "S", Amount: 6},
				{Label: "O", Amount: 24},
				{Label: "Ho", Amount: 2}},
			expected: []Atom{
				{Label: "K2O", Amount: 16.713129118300564},
				{Label: "MgO", Amount: 7.151185901417123},
				{Label: "SO3", Amount: 42.61382168343717},
				{Label: "Ho2O3", Amount: 33.52186329684514}},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			m := molarMass{tt.parsed}
			result, _ := m.oxidePercent()
			if !slices.Equal(result, tt.expected) {
				t.Errorf("molarMass() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestMolarMass_oxidePercent_noCustomOxides(t *testing.T) {
	tests := []struct {
		name     string
		parsed   []Atom
		expected []Atom
	}{
		{
			name: "no custom oxides",
			parsed: []Atom{
				{Label: "Ba", Amount: 1},
				{Label: "Fe", Amount: 1},
				{Label: "O", Amount: 4}},
			expected: []Atom{
				{Label: "BaO", Amount: 65.75731388539238},
				{Label: "Fe2O3", Amount: 34.242686114607615}},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			m := molarMass{tt.parsed}
			result, _ := m.oxidePercent()
			if !slices.Equal(result, tt.expected) {
				t.Errorf("molarMass() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestMolarMass_oxidePercent_customOxides(t *testing.T) {
	tests := []struct {
		name     string
		parsed   []Atom
		expected []Atom
	}{
		{
			name: "custom Fe oxide",
			parsed: []Atom{
				{Label: "Ba", Amount: 1},
				{Label: "Fe", Amount: 1},
				{Label: "O", Amount: 4}},
			expected: []Atom{
				{Label: "BaO", Amount: 66.51800627323722},
				{Label: "Fe3O4", Amount: 33.48199372676278}},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			m := molarMass{tt.parsed}
			result, _ := m.oxidePercent("Fe3O4")
			if !slices.Equal(result, tt.expected) {
				t.Errorf("molarMass() = %v, expected %v", result, tt.expected)
			}
		})
	}
}

func TestMolarMass_oxidePercent_wrongCustomOxides_1(t *testing.T) {
	t.Run("", func(t *testing.T) {
		data := []Atom{{Label: "Ba", Amount: 1},
			{Label: "Fe", Amount: 1},
			{Label: "O", Amount: 4}}
		m := molarMass{data}
		_, err := m.oxidePercent("Fe3O4I2")
		if err == nil {
			t.Error("want error for wrong oxide, got nil")
		}
	})
}

func TestMolarMass_oxidePercent_wrongCustomOxides_2(t *testing.T) {
	t.Run("", func(t *testing.T) {
		data := []Atom{{Label: "Ba", Amount: 1},
			{Label: "Fe", Amount: 1},
			{Label: "O", Amount: 4}}
		m := molarMass{data}
		_, err := m.oxidePercent("Fe3I2")
		if err == nil {
			t.Error("want error for wrong oxide, got nil")
		}
	})
}

/*

 */
