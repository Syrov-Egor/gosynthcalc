# gosynthcalc
<p align="center">
    <img src="data/Gopher_flask.svg" width="300" height="300">
</p>

Go library for calculating the masses of substances required for chemical synthesis directly from the reaction string. It includes solutions for all intermediate steps, including chemical formula parsing, molar mass calculation and reaction balancing with different matrix methods. This is a Go rewrite of [chemsynthcalc](https://github.com/Syrov-Egor/chemsynthcalc).

## Documentation
Detailed docs are presented for [Python version](https://syrov-egor.github.io/chemsynthcalc/). I tried to remain the structure of the project and object names are close to the Python version as possible.

## Example use
Let's say that we need to prepare 3 grams of [YBCO](https://en.wikipedia.org/wiki/Yttrium_barium_copper_oxide) by solid-state synthesis from respective carbonates. The reaction string will look something like this (to simplify, let's leave it without oxygen nonstoichiometry):
```Go
import "github.com/Syrov-Egor/gosynthcalc/"

reactionStr := "BaCO3 + Y2(CO3)3 + CuCO3 + O2 = YBa2Cu3O7 + CO2"
```
Now, we can create a chemical reaction object of the ChemicalReaction struct, which will be used in the calculation. We need to specify the arguments for our particular case:
```Go
reacOpts := ReactionOptions{
                Rmode:       Balance,
                Target:     0,
	            TargerMass: 3.0,
	            Intify:     true,
	            Precision:  8,
	            Tolerance:  1e-8,
	            }
reaction, _ := NewChemicalReaction(reactionStr, reacOpts) // Errors are supressed in this example
```
Now, to perform the automatic calculation, all we need to do is to put:
```Go
out, _ := reaction.Output()
fmt.Println(out)
```
And we get our output in the terminal:
```
	initial reaction: BaCO3+Y2(CO3)3+CuCO3+O2=YBa2Cu3O7+CO2
	reaction matrix:
	⎡1  0  0  0  2  0⎤
	⎢1  3  1  0  0  1⎥
	⎢3  9  3  2  7  2⎥
	⎢0  2  0  0  1  0⎥
	⎣0  0  1  0  3  0⎦
	mode: Balance
	formulas: [BaCO3 Y2(CO3)3 CuCO3 O2 YBa2Cu3O7 CO2]
	coefficients: [8 2 12 1 4 26]
	coefficients normalized: [2 0.5 3 0.25 1 6.5]
	algorithm: inverse
	is balanced: true
	final reaction: 8BaCO3+2Y2(CO3)3+12CuCO3+O2=4YBa2Cu3O7+26CO2
	final reaction normalized: 2BaCO3+0.5Y2(CO3)3+3CuCO3+0.25O2=YBa2Cu3O7+6.5CO2
	molar masses: [197.335 357.8357 123.554 31.998 666.1908 44.009]
	target: YBa2Cu3O7
	masses: [1.7773 0.8057 1.6692 0.036 3 1.2882]
	BaCO3      M = 197.335   g/mol  m = 1.7773  g
	Y2(CO3)3   M = 357.8357  g/mol  m = 0.8057  g
	CuCO3      M = 123.554   g/mol  m = 1.6692  g
	O2         M = 31.998    g/mol  m = 0.036   g
	YBa2Cu3O7  M = 666.1908  g/mol  m = 3       g
	CO2        M = 44.009    g/mol  m = 1.2882  g
```

## Features
* Formula parsing
```Go
form, _ := NewChemicalFormula("C2H5OH")
fmt.Println(form.ParsedFormula())
//['C': 2 'H': 6 'O': 1]
```
* Calculation of the molar mass 
```Go
form, _ := NewChemicalFormula("C2H5OH")
fmt.Println(form.MolarMass())
//46.069
```
* [Mass](https://en.wikipedia.org/wiki/Mass_fraction_(chemistry)), [atomic](https://en.wikipedia.org/wiki/Mole_fraction), and [oxide](https://d32ogoqmya1dw8.cloudfront.net/files/introgeo/studio/examples/minex02.pdf) percent calculations (including user-defined oxides).
```Go
form, _ := NewChemicalFormula("C2H5OH")
fmt.Println(form.MassPercent())
fmt.Println(form.AtomicPercent())
fmt.Println(form.OxidePercent())
//['C': 52.14352384 'H': 13.12813389 'O': 34.72834227]
//['C': 22.22222222 'H': 66.66666667 'O': 11.11111111]
//['CO2': 61.95701907 'H2O': 38.04298093]
```
* Auto-balancing chemical equations by 4 different matrix methods in `Balance` mode:
```Go
reacStr := "K4Fe(CN)6 + KMnO4 + H2SO4 = KHSO4 + Fe2(SO4)3 + MnSO4 + HNO3 + CO2 + H2O"
reac, _ := NewChemicalReaction(reacStr)
fmt.Println(reac.FinalReaction())
//10K4Fe(CN)6+122KMnO4+299H2SO4=162KHSO4+5Fe2(SO4)3+122MnSO4+60HNO3+60CO2+188H2O
```
* Calculation of masses for user-defined coefficients in `Force` (calculates regardless of balance) and `Check` (checks if reaction is balanced by user-defined coefficients) modes.
```Go
reacOpts := ReactionOptions{
		Rmode:      Force,
		Target:     0,
		TargerMass: 1.0,
		Intify:     true,
		Precision:  8,
		Tolerance:  1e-8,
	}
reaction, _ := NewChemicalReaction("BaCO3+TiO2=BaTiO3", reacOpts) //We can drop CO2 product and still get masses in this mode. 
fmt.Println(reaction.Masses())
//[0.84623763 0.34248749 1]
```
```Go
reacOpts := ReactionOptions{
		Rmode:      Check,
		Target:     0,
		TargerMass: 1.0,
		Intify:     true,
		Precision:  8,
		Tolerance:  1e-8,
	}
reaction, _ := NewChemicalReaction("H2+O2=2H2O", reacOpts) //Obviously not balanced
coefs, err := reaction.Coefficients()
if err != nil {
	fmt.Println(err)
	return
	}
fmt.Println(coefs.Result)
//Reaction is not balanced
```
* Calculation of coefficients individually by each of 4 different algorithms (inverse, general pseudoinverse, partial pseudoinverse and combinatorial algorithms).

## License
The code is provided under the MIT license.
The Go gopher was designed by [Renee French](http://reneefrench.blogspot.com/).
[gopher.svg](https://github.com/golang-samples/gopher-vector) was created by [Takuya Ueda](https://twitter.com/tenntenn).

## Contact
If you have any questions, please contact **Egor Syrov** at syrov_ev@mail.ru or
create an issue at github https://github.com/Syrov-Egor/gosynthcalc/issues.
