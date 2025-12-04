// # gosynthcalc
//
// Go library for calculating the masses of substances required for chemical synthesis directly from the reaction string.
// It includes solutions for all intermediate steps, including chemical formula parsing, molar mass calculation and reaction
// balancing with different matrix methods. This is rewrite of [Python chemsynthcalc project] in Go, so for the deep detailed documentation
// one should check [chemsynthcalc docs]. The main purpose of rewrite is to use it with [Wails] for GUIs.
//
// # Usage
//
// Let's say that we need to prepare 3 grams of [YBCO] by solid-state synthesis from respective carbonates.
// The reaction string will look something like this (to simplify, let's leave it without oxygen nonstoichiometry):
//
//	import g "github.com/Syrov-Egor/gosynthcalc/"
//
//	reactionStr := "BaCO3 + Y2(CO3)3 + CuCO3 + O2 = YBa2Cu3O7 + CO2"
//
// Now, we can create a chemical reaction object of the ChemicalReaction struct, which will be used in the calculation. We need to specify the arguments for our particular case :
//
//	reacOpts := g.ReactionOptions{
//	                Rmode:       Balance,
//	                Target:     0,
//	                TargerMass: 3.0,
//	                Intify:     true,
//	                Precision:  8,
//	                Tolerance:  1e-8,
//	            }
//
//	reaction, _ := g.NewChemicalReaction(reactionStr, reacOpts) // Errors are supressed in this example
//
// Now, to perform the automatic calculation, all we need to do is to put:
//
//	out, _ := reaction.Output()
//	fmt.Println(out)
//
// And we get our output in the terminal:
//
//	initial reaction: BaCO3+Y2(CO3)3+CuCO3+O2=YBa2Cu3O7+CO2
//	reaction matrix:
//	⎡1  0  0  0  2  0⎤
//	⎢1  3  1  0  0  1⎥
//	⎢3  9  3  2  7  2⎥
//	⎢0  2  0  0  1  0⎥
//	⎣0  0  1  0  3  0⎦
//	mode: Balance
//	formulas: [BaCO3 Y2(CO3)3 CuCO3 O2 YBa2Cu3O7 CO2]
//	coefficients: [8 2 12 1 4 26]
//	coefficients normalized: [2 0.5 3 0.25 1 6.5]
//	algorithm: inverse
//	is balanced: true
//	final reaction: 8BaCO3+2Y2(CO3)3+12CuCO3+O2=4YBa2Cu3O7+26CO2
//	final reaction normalized: 2BaCO3+0.5Y2(CO3)3+3CuCO3+0.25O2=YBa2Cu3O7+6.5CO2
//	molar masses: [197.335 357.8357 123.554 31.998 666.1908 44.009]
//	target: YBa2Cu3O7
//	masses: [1.7773 0.8057 1.6692 0.036 3 1.2882]
//	BaCO3      M = 197.335   g/mol  m = 1.7773  g
//	Y2(CO3)3   M = 357.8357  g/mol  m = 0.8057  g
//	CuCO3      M = 123.554   g/mol  m = 1.6692  g
//	O2         M = 31.998    g/mol  m = 0.036   g
//	YBa2Cu3O7  M = 666.1908  g/mol  m = 3       g
//	CO2        M = 44.009    g/mol  m = 1.2882  g
//
// [YBCO]: https://en.wikipedia.org/wiki/Yttrium_barium_copper_oxide
// [Python chemsynthcalc project]: https://github.com/Syrov-Egor/chemsynthcalc
// [chemsynthcalc docs]: https://syrov-egor.github.io/chemsynthcalc/API/
// [Wails]: https://wails.io/
package gosynthcalc

import (
	"github.com/Syrov-Egor/gosynthcalc/internal/chemformula"
	"github.com/Syrov-Egor/gosynthcalc/internal/chemreaction"
)

// A struct for operations on a single chemical formula.
// It should be constructed with [NewChemicalFormula] and can calculate
// parsed formula, molar mass, mass percent, atomic percent,
// oxide percent.
type ChemicalFormula = chemformula.ChemicalFormula

// A struct for operations on a single chemical formula.
// It should be constructed with [NewChemicalReaction] and can calculate
// coefficients of reaction and output masses of compounds.
type ChemicalReaction = chemreaction.ChemicalReaction

// Defaults are
//
//	ReacOptions{
//					Rmode:       Balance,
//					Target:     0,
//					TargerMass: 1.0,
//					Intify:     true,
//					Precision:  8,
//					Tolerance:  1e-8,
//				}
//
//	- Rmode: Coefficients calculation mode [ReactionMode]
//	- Target: Index of target compound (0 by default, or first compound in the products), can be negative (limited by reactant)
//	- Target_mass: Desired mass of target compound (in grams)
//	- Intify: Is it required to convert the coefficients to integer values?
//	- Precision: Value of rounding precision (8 by default)
//	- Tolerance: Tolerance for comparing floats (1e-8 by default)
type ReactionOptions = chemreaction.ReacOptions

//  1. The "force" mode is used when a user enters coefficients
//     in the reaction string and wants the masses to be calculated
//     whether the reaction is balanced or not.
//
//  2. "check" mode is the same as force, but with reaction
//     balance checks.
//
//  3. "balance" mode  tries to automatically calculate
//     coefficients from the reaction string.
type ReactionMode = chemreaction.Mode

const (
	Force   ReactionMode = chemreaction.Force
	Check   ReactionMode = chemreaction.Check
	Balance ReactionMode = chemreaction.Balance
)

type MethodResult = chemreaction.MethodResult

// Builder function to create [ChemicalFormula] object.
func NewChemicalFormula(formula string, precision ...uint) (*ChemicalFormula, error) {
	return chemformula.NewChemicalFormula(formula, precision...)
}

// Builder function to create [ChemicalReaction] object.
func NewChemicalReaction(reaction string, options ...ReactionOptions) (*ChemicalReaction, error) {
	return chemreaction.NewChemicalReaction(reaction, options...)
}
