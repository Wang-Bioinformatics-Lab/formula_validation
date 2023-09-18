#!/usr/bin/env python

"""
This Python module contains the tests of the class Formula. 

@contents :  This Python module contains the tests of the class Formula. 
this Formula class. It uses all the chemcalc features to process complex molecular formulas as the one with nested parenthesis and its possibility of checking the possible molecular formulas from a monoisotopic mass taking into account the number of elements
@project :  N/A
@program :  N/A
@file :  Formula.py
@author :  Alberto Gil De la Fuente (alberto.gilf@gmail.com)

@version :  0.0.1, 20 July 2023
@information :  The Zen of Python
          https://www.python.org/dev/peps/pep-0020/
        Style Guide for Python Code
          https://www.python.org/dev/peps/pep-0008/
        Example NumPy Style Python Docstrings
          http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html
        doctest – Testing through documentation
          https://pymotw.com/2/doctest/

@copyright :  Copyright 2023 GNU AFFERO GENERAL PUBLIC.
        All rights are reserved. Reproduction in whole or in part is
        prohibited without the written consent of the copyright owner.
"""

__author__    = "Alberto Gil de la Fuente"
__copyright__   = "GPL License version 3"

import unittest
from formula_validation.Formula import Formula
from formula_validation.Element import Element_type, element_weights
from formula_validation.IncorrectFormula import IncorrectFormula
from formula_validation.NotFoundElement import NotFoundElement


class TestFormula(unittest.TestCase):
  def setUp(self):
    # Set up any common resources needed for tests
    pass

  def tearDown(self):
    # Clean up after each test
    pass

  def test_constructor_from_dictionary(self):
    my_elements = {'C': 5, 'H': 5, 'N': 4}
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula(my_elements, adduct)
    self.assertTrue(True)

  def test_monoisotopic_mass(self):
    my_elements = {'C': 5, 'H': 5, 'N': 4}
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula(my_elements, adduct)
    expected_monoisotopic_mass = 121.05142
    formula_monoisotopic_mass = my_formula.get_monoisotopic_mass()
    self.assertAlmostEqual(expected_monoisotopic_mass, formula_monoisotopic_mass, delta=1e-3)

  def test_formula_with_adduct(self):
    my_elements = {Element_type.C: 5, Element_type.H: 5, Element_type.N: 4}
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula(my_elements, adduct)
    self.assertTrue(True)

  def test_monoisotopic_mass_with_adduct(self):
    my_elements = {Element_type.C: 5, Element_type.H: 5, Element_type.N: 4}
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula(my_elements, adduct)
    expected_monoisotopic_mass = 122.05862
    formula_monoisotopic_mass = my_formula.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(expected_monoisotopic_mass, formula_monoisotopic_mass, delta=1e-4)

  def test_formula_from_str_plain(self):
    my_formula_str = "H5C5N4ONaK"
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_str(my_formula_str, adduct)
    elements_expected = {Element_type.H: 5, Element_type.C: 5, Element_type.N: 4, Element_type.O: 1, Element_type.Na: 1, Element_type.K: 1}
    self.assertEqual(my_formula.get_elements(),elements_expected)
  
  def test_formula_from_str_with_parenthesis(self):
    my_formula_str = "H5C5N4O(H2O)2NaK"
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_str(my_formula_str, adduct)
    elements_expected = {Element_type.H: 9, Element_type.C: 5, Element_type.N: 4, Element_type.O: 3, Element_type.Na: 1, Element_type.K: 1}
    self.assertEqual(my_formula.get_elements(),elements_expected)
    
  def test_formula_from_smiles_1(self):
    smiles = 'CC(C)CC1NC(=O)C(C)NC(=O)C(=C)N(C)C(=O)CCC(NC(=O)C(C)C(NC(=O)C(CCCNC(N)=N)NC(=O)C(C)C(NC1=O)C(O)=O)\\C=C\\C(\\C)=C\\C(C)C(O)Cc1ccccc1)C(O)=O'
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_smiles(smiles, adduct)
    elements_expected = {Element_type.H: 72, Element_type.C: 48, Element_type.N: 10, Element_type.O: 12}
    self.assertEqual(my_formula.get_elements(),elements_expected)

  def test_formula_from_smiles_2(self):
    smiles = 'CCC[C@@H](C)[C@@H]([C@H](C)[C@@H]1[C@H]([C@H](Cc2nc(cs2)C3=N[C@](CS3)(C4=N[C@](CS4)(C(=O)N[C@H]([C@H]([C@H](C(=O)O[C@H](C(=O)N[C@H](C(=O)O1)[C@@H](C)O)[C@@H](C)CC)C)O)[C@@H](C)CC)C)C)OC)C)O'
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_smiles(smiles, adduct)
    elements_expected = {Element_type.H: 73, Element_type.C: 45, Element_type.N: 5, Element_type.O: 10, Element_type.S: 3}
    self.assertEqual(my_formula.get_elements(),elements_expected)

  def test_formula_from_smiles_3(self):
    smiles = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_smiles(smiles, adduct)
    elements_expected = {Element_type.H: 38, Element_type.C: 24, Element_type.N: 1, Element_type.O: 5, Element_type.Cl: 1}
    self.assertEqual(my_formula.get_elements(),elements_expected)

  def test_monoisotopic_mass(self):
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = 'None'
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    expected_value = 0.0228
    current_value = Formula.ppm_to_absolute(my_formula.get_monoisotopic_mass())
    self.assertAlmostEqual(current_value, expected_value, delta=1e-3)

  def test_formula_and_external_molecular_mass(self):
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    expected_value = True
    external_mass=455.24
    current_value = my_formula.check_monoisotopic_mass(external_mass,50)
    self.assertAlmostEqual(current_value, expected_value, delta=1e-3)


  def test_monoisotopic_mass_with_adduct(self):
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[M+C2H2O-H]-' 
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    expected_value =496.24713858
    current_value = my_formula.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(current_value, expected_value, delta=1e-3)

  def test_monoisotopic_mass_with_adduct_2(self):
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[M-2H2O+2H]2+' 
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    expected_value =210.61863642
    current_value = my_formula.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(current_value, expected_value, delta=1e-3)


  def test_monoisotopic_mass_with_adduct_3(self):
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[5M+Ca]2+' 
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    expected_value = 1158.09037142
    current_value = my_formula.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(current_value, expected_value, delta=1e-3)

  def test_ppm_differences_with_adduct(self):
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[5M+Ca]2+' 
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    experimental_mass = 1158.099
    expected_value = 7.45
    current_value = my_formula.ppm_difference_with_exp_mass(experimental_mass)
    self.assertAlmostEqual(current_value, expected_value, delta=1e-3)

  def test_addition_two_formulas(self):
    formula_1 = 'C4H4O2'
    formula_2 = 'C4H6O1'
    adduct = None 
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_2 = Formula.formula_from_str(formula_2, adduct)
    expected_value = Formula.formula_from_str('C8H10O3', adduct)
    current_value = my_formula_1 + my_formula_2
    self.assertEqual(current_value, expected_value)

  def test_subtraction_two_formulas(self):
    formula_1 = 'C4H4O2'
    formula_2 = 'C4H2O1'
    adduct = None 
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_2 = Formula.formula_from_str(formula_2, adduct)
    expected_value = Formula.formula_from_str('H2O', adduct)
    current_value = my_formula_1 - my_formula_2
    self.assertEqual(current_value, expected_value)


  def test_subtraction_two_formulas_raise_error(self):
    formula_1 = 'C4H4O2'
    formula_2 = 'C4H2O1'
    adduct = None 
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_2 = Formula.formula_from_str(formula_2, adduct)
    with self.assertRaises(IncorrectFormula) as context:
      current_value = my_formula_2 - my_formula_1
        
    self.assertIn("The subtraction of these two formulas contain a negative number of ",str(context.exception))
    
  def test_verify_fragment_mz_no_adduct(self):
    formula_1 = 'C9H9O2S2'
    adduct = None 
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    current_value = my_formula_1.check_possible_fragment_mz(fragment_mz=200,ppm=50)
    expected_value = True
    self.assertEqual(current_value, expected_value)
  
  def test_verify_fragment_mz_with_adduct(self):
    formula_1 = 'C9H9O2S2'
    adduct = '[M-C2+H]+' 
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    current_value = my_formula_1.check_possible_fragment_mz(fragment_mz=200,ppm=50)
    expected_value = False
    self.assertEqual(current_value, expected_value)


  def test_percentage_intensity_fragments_explained_by_formula(self):
    
    formula_1 = 'C9H9O2S2'
    adduct = '[M+H]+' 
    import pickle
    spectra_1 = pickle.load(open("example_spectra_for_Alberto.pkl", 'rb'))
    mzs = spectra_1['m/z array']
    intensities = spectra_1['intensity array']
    fragments_mz_intensities = dict(zip(mzs, intensities))
    metadata = pickle.load(open("example_metadata_for_Alberto.pkl", 'rb'))
    smiles = metadata[1]["Smiles"]
    formula = Formula.formula_from_smiles(smiles, adduct)

    
    current_value = formula.percentage_intensity_fragments_explained_by_formula(fragments_mz_intensities=fragments_mz_intensities,ppm=50)
    expected_value = 0.886
    
    self.assertAlmostEqual(current_value, expected_value, delta=1e-2)

  def test_formula_str(self):
    expected_value = "C12H3N3O+[M+H-C2]+"
    formula_1 = 'C12H3N3O'
    adduct = '[M-C2+H]+'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_1_str = str(my_formula_1)
    self.assertEqual(my_formula_1_str,expected_value)

  def test_final_formula_with_adduct(self):
    expected_value = "[C10H4N3O]+"
    formula_1 = 'C12H3N3O'
    adduct = '[M-C2+H]+'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_1_str = my_formula_1.get_final_formula_with_adduct()
    self.assertEqual(my_formula_1_str,expected_value)

if __name__ == "__main__":
  unittest.main()