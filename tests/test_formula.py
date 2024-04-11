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

  def test_monoisotopic_mass_with_isotopes(self):
    my_elements = '[13]C2[14]C5C5[2]H6H8N4O4S'
    adduct = 'None'
    my_formula = Formula.formula_from_str(my_elements, adduct)
    expected_monoisotopic_mass = 328.133317657
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

  def test_formula_from_inchi_1(self):
    inchi = 'InChI=1S/C48H72N10O12/c1-25(2)22-36-45(66)57-39(47(69)70)29(6)41(62)54-34(16-13-21-51-48(49)50)44(65)53-33(18-17-26(3)23-27(4)37(59)24-32-14-11-10-12-15-32)28(5)40(61)55-35(46(67)68)19-20-38(60)58(9)31(8)43(64)52-30(7)42(63)56-36/h10-12,14-15,17-18,23,25,27-30,33-37,39,59H,8,13,16,19-22,24H2,1-7,9H3,(H,52,64)(H,53,65)(H,54,62)(H,55,61)(H,56,63)(H,57,66)(H,67,68)(H,69,70)(H4,49,50,51)/b18-17+,26-23+'
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_inchi(inchi, adduct)
    elements_expected = {Element_type.H: 72, Element_type.C: 48, Element_type.N: 10, Element_type.O: 12}
    self.assertEqual(my_formula.get_elements(),elements_expected)

  def test_formula_from_inchi_2(self):
    inchi = 'InChI=1S/C45H73N5O10S3/c1-14-17-24(6)34(52)26(8)37-25(7)30(58-13)18-31-46-29(19-61-31)39-49-45(12,21-62-39)43-50-44(11,20-63-43)42(57)48-32(22(4)15-2)35(53)27(9)40(55)59-36(23(5)16-3)38(54)47-33(28(10)51)41(56)60-37/h19,22-28,30,32-37,51-53H,14-18,20-21H2,1-13H3,(H,47,54)(H,48,57)/t22-,23-,24+,25-,26-,27+,28+,30-,32-,33-,34-,35-,36-,37?,44+,45?/m0/s1'
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_inchi(inchi, adduct)
    elements_expected = {Element_type.H: 73, Element_type.C: 45, Element_type.N: 5, Element_type.O: 10, Element_type.S: 3}
    self.assertEqual(my_formula.get_elements(),elements_expected)

  def test_formula_from_inchi_3(self):
    inchi = 'InChI=1S/C24H38ClNO5/c1-3-4-5-6-8-11-19(30-2)12-9-7-10-13-22(29)26-17-18(16-25)24-21(28)15-14-20(27)23(24)31-24/h7,9,16,19-20,23,27H,3-6,8,10-15,17H2,1-2H3,(H,26,29)/b9-7+,18-16-/t19-,20-,23-,24+/m0/s1'
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_inchi(inchi, adduct)
    elements_expected = {Element_type.H: 38, Element_type.C: 24, Element_type.N: 1, Element_type.O: 5, Element_type.Cl: 1}
    self.assertEqual(my_formula.get_elements(),elements_expected)

  def test_formula_from_inchi_4(self):
    inchi = 'InChI=1S/C3H10N2/c4-2-1-3-5/h1-5H2'
    my_formula = Formula.formula_from_inchi(inchi)
    elements_expected = {Element_type.H: 10, Element_type.C: 3, Element_type.N: 2}
    self.assertEqual(my_formula.get_elements(),elements_expected)
  
  def test_formula_from_inchi_with_a_charge(self):
    inchi = 'InChI=1S/C5H14NO/c1-6(2,3)4-5-7/h7H,4-5H2,1-3H3/q+1'
    my_formula = Formula.formula_from_inchi(inchi)
    elements_expected = {Element_type.H: 14, Element_type.C: 5, Element_type.N: 1, Element_type.O: 1}
    self.assertEqual(my_formula.get_elements(),elements_expected)

  def test_formula_from_inchi_with_a_proton(self):
    inchi = 'InChI=1S/C26H24O14/c1-35-14-3-9(4-15(36-2)19(14)30)23-24(40-26-22(33)21(32)20(31)16(8-28)39-26)17-11(7-27)25(34)38-13-6-10(29)5-12(37-23)18(13)17/h3-7,16,20-22,26,28,30-34H,8H2,1-2H3/p+1'
    my_formula = Formula.formula_from_inchi(inchi)
    elements_expected = {Element_type.H: 14, Element_type.C: 5, Element_type.N: 1, Element_type.O: 1}
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

  def test_final_formula_with_no_adduct(self):
    expected_value = "C5H4O4"
    formula_1 = 'C5H4O4'
    adduct = None
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_1_str = my_formula_1.get_final_formula_with_adduct()
    self.assertEqual(my_formula_1_str,expected_value)

  def test_final_formula_with_adduct_neutral(self):
    expected_value = "C5H4O4Na"
    formula_1 = 'C5H5O4'
    adduct = '[M-H+Na]'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_1_str = my_formula_1.get_final_formula_with_adduct()
    self.assertEqual(my_formula_1_str,expected_value)

  def test_final_formula_with_multimer2(self):
    expected_value = "C10H9O8Na"
    formula_1 = 'C5H5O4'
    adduct = '[2M-H+Na]'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_1_str = my_formula_1.get_final_formula_with_adduct()
    self.assertEqual(my_formula_1_str,expected_value)
  
  
  def test_natively_pos_charged_molecule(self):
    expected_value = "C10H9O8+"
    formula_1 = 'C10H9O8+'
    adduct = 'None'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_1_str = my_formula_1.get_final_formula_with_adduct()
    self.assertEqual(my_formula_1_str,expected_value)

  def test_natively_pos_charged_molecule_mass(self):
    expected_value = 257.02919
    formula_1 = 'C10H9O8+'
    adduct = 'None'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    mz = my_formula_1.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(mz,expected_value, delta=1e-3)


  def test_natively_neg_charged_molecule(self):
    expected_value = "C10H9O8-"
    formula_1 = 'C10H9O8-'
    adduct = 'None'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_1_str = my_formula_1.get_final_formula_with_adduct()
    self.assertEqual(my_formula_1_str,expected_value)
  
  def test_natively_neg_charged_molecule_mass(self):
    expected_value = 257.03029
    formula_1 = 'C10H9O8-'
    adduct = 'None'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    mz = my_formula_1.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(mz,expected_value, delta=1e-3)
  
  def test_natively_pos_charged_molecule_with_neutral_adduct(self):
    expected_value = "C20H17O16Na+"
    formula_1 = 'C10H9O8+'
    adduct = '[2M-H+Na]'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_1_str = my_formula_1.get_final_formula_with_adduct()
    self.assertEqual(my_formula_1_str,expected_value)

  def test_natively_pos_charged_molecule_with_neutral_adduct_mass(self):
    expected_value = 536.04088
    formula_1 = 'C10H9O8+'
    adduct = '[2M-H+Na]'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    mz = my_formula_1.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(mz,expected_value, delta=1e-3)

  def test_natively_pos_charged_molecule_with_charged_adduct(self):
    expected_value = "[C20H18O16Na]+2"
    formula_1 = 'C10H9O8+'
    adduct = '[2M+Na]+'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_1_str = my_formula_1.get_final_formula_with_adduct()
    self.assertEqual(my_formula_1_str,expected_value)
  
  def test_natively_pos_charged_molecule_with_charged_adduct_mass(self):
    expected_value = 268.52408
    formula_1 = 'C10H9O8+'
    adduct = '[2M+Na]+'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    mz = my_formula_1.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(mz,expected_value, delta=1e-3)

  def test_natively_pos_charged_molecule_with_charged_adduct_with_parenthesis(self):
    expected_value = "[C20H18O16Na]+2"
    formula_1 = 'C10H9O8(+)'
    adduct = '[2M+Na]+'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_1_str = my_formula_1.get_final_formula_with_adduct()
    self.assertEqual(my_formula_1_str,expected_value)

  def test_natively_pos_charged_molecule_with_charged_adduct_mass_with_parenthesis(self):
    expected_value = 268.52408
    formula_1 = 'C10H9O8(+1)'
    adduct = '[2M+Na]+'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    mz = my_formula_1.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(mz,expected_value, delta=1e-3)

  def test_natively_pos_charged_molecule_adduct_M_plus(self):
    expected_value = 257.029197
    formula_1 = 'C10H9O8(+1)'
    adduct = '[M]+'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    mz = my_formula_1.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(mz,expected_value, delta=1e-3)
  
  def test_natively_pos_charged_molecule_adduct_M(self):
    expected_value = 257.029197
    formula_1 = 'C10H9O8(+1)'
    adduct = '[M]'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    mz = my_formula_1.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(mz,expected_value, delta=1e-3)

  def test_natively_pos_charged_molecule_adduct_M_minus(self):
    expected_value = 257.029197
    formula_1 = 'C10H9O8(+1)'
    adduct = '[M]'
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    mz = my_formula_1.get_monoisotopic_mass_with_adduct()
    self.assertAlmostEqual(mz,expected_value, delta=1e-3)
  
if __name__ == "__main__":
  unittest.main()