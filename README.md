# Formula Validation Python Module

This Python module contains a `Formula` class for working with chemical formulas, as well as methods for creating, manipulating, and analyzing formulas. It also includes functionality for dealing with adducts and calculating monoisotopic masses.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
  - [Creating Formula Objects](#creating-formula-objects)
  - [Basic Operations](#basic-operations)
  - [Calculating Mass](#calculating-mass)
  - [Fragment Analysis](#fragment-analysis)
- [Examples](#examples)
- [License](#license)

## Introduction

This Python module provides a `Formula` class that allows you to work with chemical formulas. It includes the following features:

- Create Formula objects from Hill notation, SMILES, and InChI.
- Perform basic mathematical operations on formulas (addition, subtraction, multiplication).
- Calculate the monoisotopic mass of a formula.
- Check if a given mass is within a specified tolerance of the formula's mass.
- Analyze possible fragment masses explained by a formula and adduct.

## Installation

To use this module, you'll need Python 3.x and the required dependencies. You can install the dependencies using pip:

```bash
pip install rdkit urllib3 rpy2
```
You should have R installed with the package devtools and 
```R
install.packages("devtools")
devtools::install_github("mjhelf/MassTools")

library(MassTools)
```
## Usage

### Creating Formula Objects

You can create a Formula object using various methods:

- `Formula.formula_from_str_hill(formula_str: str, adduct: str) -> 'Formula'`: Create a Formula object from a chemical formula string in Hill notation.

- `Formula.formula_from_str(formula_str: str, adduct: str, no_api: bool = False) -> 'Formula'`: Create a Formula object from a chemical formula string. You can disable API calls for formula resolution by setting `no_api` to `True`.

- `Formula.formula_from_smiles(smiles: str, adduct: str, no_api: bool = False) -> 'Formula'`: Create a Formula object from a SMILES string representing a molecular structure.

- `Formula.formula_from_inchi(inchi: str, adduct: str, no_api: bool = False) -> 'Formula'`: Create a Formula object from an InChI string representing a molecular structure.

### Basic Operations

You can perform various operations on Formula objects:

- Addition: `formula1 + formula2`
- Subtraction: `formula1 - formula2`
- Multiplication: `formula * num`

### Calculating Mass

You can calculate the monoisotopic mass of a formula and check if it matches an external mass:

- `get_monoisotopic_mass() -> float`: Get the monoisotopic mass of the formula.
- `get_monoisotopic_mass_with_adduct() -> float`: Get the monoisotopic mass of the formula, considering the adduct.
- `check_monoisotopic_mass(external_mass: Union[float, int], mass_tolerance_in_ppm: Union[int, float] = __default_ppm) -> bool`: Check if the monoisotopic mass is within a specified tolerance of an external mass.
- `check_monoisotopic_mass_with_adduct(external_mass: Union[float, int], mass_tolerance_in_ppm: Union[int, float] = __default_ppm) -> bool`: Check if the monoisotopic mass, considering the adduct, is within a specified tolerance of an external mass.

### Fragment Analysis

You can analyze potential fragment masses explained by a formula and adduct:

- `check_possible_fragment_mz(fragment_mz: Union[float, int], ppm: Union[float, int] = __default_ppm) -> bool`: Check if a fragment mass can be explained by the formula and adduct.

- `percentage_intensity_fragments_explained_by_formula(fragments_mz_intensities: Dict[Union[float, int], Union[float, int]], ppm: Union[float, int] = __default_ppm) -> float`: Calculate the percentage of intensity of fragments explained by the formula and adduct.

## Examples

Here are some examples of how to use the Formula class:

```python
# Create Formula objects
formula1 = Formula.formula_from_str_hill("C5H5O4", "[M+H]+")
formula2 = Formula.formula_from_smiles("CCO", "[M+NH4]+")
formula3 = formula1 + formula2

# Calculate monoisotopic mass
mass1 = formula1.get_monoisotopic_mass()
mass2 = formula2.get_monoisotopic_mass_with_adduct()
print(f"Mass of formula1: {mass1}")
print(f"Mass of formula2 with adduct: {mass2}")

# Check mass against an external mass with a tolerance of 5 ppm
check_monoisotopic_mass = formula1.check_monoisotopic_mass(121.05142,5)
check_monoisotopic_mass_with_adduct = formula1.check_monoisotopic_mass_with_adduct(122.05862,5)
```
