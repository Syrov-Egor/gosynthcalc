package formula

type MolarMass struct {
	parsed []Atom
}

func (m MolarMass) atomicMasses() []float64 {
	masses := make([]float64, len(m.parsed))
	for i, atom := range m.parsed {
		masses[i] = PeriodicTable[atom.Label].Weight * atom.Amount
	}
	return masses
}

func (m MolarMass) calcMolarMass() float64 {
	var sum float64 = 0.0
	for _, mass := range m.atomicMasses() {
		sum += mass
	}
	return sum
}
