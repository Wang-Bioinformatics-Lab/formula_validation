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
from Formula import Formula
from Element import Element_type, element_weights
from IncorrectFormula import IncorrectFormula
from NotFoundElement import NotFoundElement


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



def main():
  
  import math

  adduct = None
  print("=================================================================.")
  print("Test Case 1A: Creating a formula from a list of string elements")
  print("=================================================================.")
  try:
    my_elements = {'C': 5, 'H': 5, 'N': 4}
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula(my_elements, adduct)
    print("Test PASSED. The method constructor from a collection of elements has been implemented correctly.")
    
    print("=================================================================.")
    print("Test Case 2A: Calculating the monoisotopic mass of a formula")
    print("=================================================================.")
    expected_monoisotopic_mass = 121.05142
    formula_monoisotopic_mass = my_formula.get_monoisotopic_mass()
    if math.isclose(formula_monoisotopic_mass,expected_monoisotopic_mass,abs_tol=0.0001):
      print("Test PASSED. The method Formula.get_monoisotopic_mass() has been implemented correctly.")
    else:
      print("Test FAILED. Check the method Formula.get_monoisotopic_mass(). Expected: " + str(expected_monoisotopic_mass) + " ACTUAL: " + str(formula_monoisotopic_mass))
  except IncorrectFormula as incfor:
      print("Test FAILED. Check the constructor of formula " + str(incfor))
  except NotFoundElement as nfe:
      print("Test FAILED. Check the map of the elements. Something is not an element " + str(nfe))

  print("=================================================================.")
  print("Test Case 1B: Creating a formula from a list of element types")
  print("=================================================================.")
  try:
    my_elements = {Element_type.C: 5, Element_type.H: 5, Element_type.N: 4}
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula(my_elements, adduct)
    print("Test PASSED. The method constructor from a collection of elements has been implemented correctly.")
    
    print("=================================================================.")
    print("Test Case 2B: Calculating the monoisotopic mass of a formula")
    print("=================================================================.")
    expected_monoisotopic_mass = 121.05142
    formula_monoisotopic_mass = my_formula.get_monoisotopic_mass()
    if math.isclose(formula_monoisotopic_mass,expected_monoisotopic_mass,abs_tol=0.0001):
      print("Test PASSED. The method Formula.get_monoisotopic_mass() has been implemented correctly.")
    else:
      print("Test FAILED. Check the method Formula.get_monoisotopic_mass(). Expected: " + str(expected_monoisotopic_mass) + " ACTUAL: " + str(formula_monoisotopic_mass))
  except IncorrectFormula as incfor:
      print("Test FAILED. Check the constructor of formula " + str(incfor))
  except NotFoundElement as nfe:
      print("Test FAILED. Check the map of the elements. Something is not an element " + str(nfe))



  print("=================================================================.")
  print("Test Case 3A: Creating a formula from a string")
  print("=================================================================.")
  try:
    my_formula_str = "H5C5N4ONaK"
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_str(my_formula_str, adduct)
    elements_expected = {Element_type.H: 5, Element_type.C: 5, Element_type.N: 4, Element_type.O: 1, Element_type.Na: 1, Element_type.K: 1}
    if my_formula.get_elements() == elements_expected:
      print("Test PASSED. The function to construct a formula representing a string has been implemented correctly.")
    else:
      print("Test FAILED. Check the constructor of formula from string")
      print(my_formula)
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check the constructor of formula from a string- Check the exception information. ")

  print("=================================================================.")
  print("Test Case 3B: Creating a formula from a string with parenthesis")
  print("=================================================================.")
  try:
    my_formula_str = "H5C5N4O(H2O)2NaK"
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_str(my_formula_str, adduct)
    elements_expected = {Element_type.H: 9, Element_type.C: 5, Element_type.N: 4, Element_type.O: 3, Element_type.Na: 1, Element_type.K: 1}
    if my_formula.get_elements() == elements_expected:
      print("Test PASSED. The function to construct a formula representing a string has been implemented correctly.")
    else:
      print("Test FAILED. Check the constructor of formula from string")
      print(my_formula)
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check the constructor of formula from a string. Check the exception information. ")

  
  print("=================================================================.")
  print("Test Case 4: Creating a formula from a SMILES")
  print("=================================================================.")
  try:
    smiles_1 = 'CC(C)CC1NC(=O)C(C)NC(=O)C(=C)N(C)C(=O)CCC(NC(=O)C(C)C(NC(=O)C(CCCNC(N)=N)NC(=O)C(C)C(NC1=O)C(O)=O)\\C=C\\C(\\C)=C\\C(C)C(O)Cc1ccccc1)C(O)=O'
    smiles_2 = 'CCC[C@@H](C)[C@@H]([C@H](C)[C@@H]1[C@H]([C@H](Cc2nc(cs2)C3=N[C@](CS3)(C4=N[C@](CS4)(C(=O)N[C@H]([C@H]([C@H](C(=O)O[C@H](C(=O)N[C@H](C(=O)O1)[C@@H](C)O)[C@@H](C)CC)C)O)[C@@H](C)CC)C)C)OC)C)O'
    smiles_3 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    Formula.formula_from_smiles('C/C1=C\CCC2=C[C@@H](OC2=O)[C@H](C(C)C)CC[C@]3(C)[C@H](O3)CC1', '[M-H2O+H]+')
    my_formula_2 = Formula.formula_from_smiles(smiles_2, adduct)
    my_formula_3 = Formula.formula_from_smiles(smiles_3, adduct)
    elements_expected_1 = {Element_type.H: 72, Element_type.C: 48, Element_type.N: 10, Element_type.O: 12}
    elements_expected_2 = {Element_type.H: 73, Element_type.C: 45, Element_type.N: 5, Element_type.O: 10, Element_type.S: 3}
    elements_expected_3 = {Element_type.H: 38, Element_type.C: 24, Element_type.N: 1, Element_type.O: 5, Element_type.Cl: 1}
    if my_formula.get_elements() == elements_expected_1:
      print("Test PASSED. The function to construct a formula by a SMILES has been implemented correctly.")
    else:
      print("Test FAILED. Check the constructor of formula from string: ")
      print("EXPECTED" + my_formula)
      print("ACTUAL" + my_formula.get_elements())
      print(my_formula)
    if my_formula_2.get_elements() == elements_expected_2:
      print("Test PASSED. The function to construct a formula by a SMILES has been implemented correctly.")
    else:
      print("Test FAILED. Check the constructor of formula from string: ")
      print("EXPECTED" + my_formula)
      print("ACTUAL" + my_formula.get_elements())
      print(my_formula)
    if my_formula_3.get_elements() == elements_expected_3:
      print("Test PASSED. The function to construct a formula by a SMILES has been implemented correctly.")
    else:
      print("Test FAILED. Check the constructor of formula from string: ")
      print("EXPECTED" + my_formula)
      print("ACTUAL" + my_formula.get_elements())
      print(my_formula)
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check the constructor of formula from a string. Check the exception information. ")

  print("=================================================================.")
  print("Test Case 5: Calculating absolute values from ppm")
  print("=================================================================.")
  try:
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = 'None'
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    expected_value = 0.0228
    current_value = Formula.ppm_to_absolute(my_formula.get_monoisotopic_mass())
    if(math.isclose(current_value,expected_value,abs_tol=0.001)):
      print("Test PASSED. The function to calculate the absolute values from ppms is correct.")
    else:
      print("Test FAILED. Check the function to calculate absolute values from ppms")
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check the function to calculate absolute values from ppms. Check the exception information. ")

  print("=================================================================.")
  print("Test Case 6: Check if a mass from a formula is less than 50 ppm ")
  print("=================================================================.")
  try:
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[M+C2H2O-H]-'
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    expected_value = True
    external_mass=455.24
    current_value = my_formula.check_monoisotopic_mass(external_mass,50)
    if(math.isclose(current_value,expected_value,abs_tol=0.001)):
      print("Test PASSED. The function to check if a formula and an external mass have a lower difference than the ppm establised of 50.")
    else:
      print("Test FAILED. Check the function to check if a formula and an external mass have a lower difference than the ppm establised. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check the function to check a formula and an external mass have a lower difference than the ppm establised. Check the exception information. ")

  print("=================================================================.")
  print("Test Case 7: Calculating the monoisotopic mass taking into account the adducts. Simple charged negative")
  print("=================================================================.")
  try:
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[M+C2H2O-H]-' 
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    expected_value =496.24713858
    current_value = my_formula.get_monoisotopic_mass_with_adduct()
    if(math.isclose(current_value,expected_value,abs_tol=0.001)):
      print("Test PASSED. The function to calculate the monoisotopic mass with an adduct in a single negative charge.")
    else:
      print("Test FAILED. Check the function to calculate the monoisotopic mass with an adduct in a single negative charge. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check the function to calculate the monoisotopic mass with an adduct in a single negative charge. Check the exception information. ")
  # '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+' 
  
  print("=================================================================.")
  print("Test Case 8: Calculating the monoisotopic mass taking into account the adducts. Double charged [M-2H2O+2H]2+")
  print("=================================================================.")
  try:
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[M-2H2O+2H]2+' 
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    expected_value =210.61863642
    current_value = my_formula.get_monoisotopic_mass_with_adduct()
    if(math.isclose(current_value,expected_value,abs_tol=0.001)):
      print("Test PASSED. The function to calculate the monoisotopic mass with an adduct in a double positive charge.")
    else:
      print("Test FAILED. Check the function to calculate the monoisotopic mass with an adduct in a double positive charge. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check the function to calculate the monoisotopic mass with an adduct in a single positive charge. Check the exception information. ")

  print("=================================================================.")
  print("Test Case 9: Calculating the monoisotopic mass taking into account the adducts. Double charged [5M+Ca]2+ and multimer 5")
  print("=================================================================.")
  try:
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[5M+Ca]2+' 
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    expected_value = 1158.09037142
    current_value = my_formula.get_monoisotopic_mass_with_adduct()
    if(math.isclose(current_value,expected_value,abs_tol=0.001)):
      print("Test PASSED. The function to calculate the monoisotopic mass with an adduct in a double positive charge and multimer 5.")
    else:
      print("Test FAILED. Check the function to calculate the monoisotopic mass with an adduct in a double positive charge and multimer 5. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check the function to calculate the monoisotopic mass with an adduct in a single positive charge and multimer 5. Check the exception information. ")

  print("=================================================================.")
  print("Test Case 10: Calculating the ppm difference between a formula and an experimental value")
  print("=================================================================.")
  try:
    smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
    adduct = '[5M+Ca]2+' 
    my_formula = Formula.formula_from_smiles(smiles_1, adduct)
    experimental_mass = 1158.099
    expected_value = 7.45
    current_value = my_formula.ppm_difference_with_exp_mass(experimental_mass)
    if(math.isclose(current_value,expected_value,abs_tol=0.01)):
      print("Test PASSED. The function Calculating the ppm difference between a formula and an experimental value.")
    else:
      print("Test FAILED. Check the function Calculating the ppm difference between a formula and an experimental value. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check the function Calculating the ppm difference between a formula and an experimental value. Check the exception information. ")

  print("=================================================================.")
  print("Test Case 11: Testing the addition of two formulas")
  print("=================================================================.")
  try:
    formula_1 = 'C4H4O2'
    formula_2 = 'C4H6O1'
    adduct = None 
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_2 = Formula.formula_from_str(formula_2, adduct)
    expected_value = Formula.formula_from_str('C8H10O3', adduct)
    current_value = my_formula_1 + my_formula_2
    if(current_value == expected_value):
      print("Test PASSED. The method to add two formulas is working.")
    else:
      print("Test FAILED. Check the method The function to add two formulas. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check the method The function to add two formulas. Check the exception information. ")

  print("=================================================================.")
  print("Test Case 12: Testing the subtraction of two formulas with a correct formula")
  print("=================================================================.")
  try:
    formula_1 = 'C4H4O2'
    formula_2 = 'C4H2O1'
    adduct = None 
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_2 = Formula.formula_from_str(formula_2, adduct)
    expected_value = Formula.formula_from_str('H2O', adduct)
    current_value = my_formula_1 - my_formula_2
    if(current_value == expected_value):
      print("Test PASSED. The method to subtract two formulas correct is working.")
    else:
      print("Test FAILED. Check The method to subtract two formulas with proper result. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. The method to subtract two formulas with proper result. Check the exception information. ")


  print("=================================================================.")
  print("Test Case 13: Testing the subtraction of two formulas with netagive elements as a result")
  print("=================================================================.")
  try:
    formula_1 = 'C4H4O2'
    formula_2 = 'C4H6O1'
    adduct = None 
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    my_formula_2 = Formula.formula_from_str(formula_2, adduct)
    expected_value = Formula.formula_from_str('C8H10O3', adduct)
    current_value = my_formula_1 - my_formula_2
    if(current_value == expected_value):
      print("Test FAILED. Check The method to subtract two formulas with an exception because the second one contains more elements of a specific type than the first one. Expected an exception, but getting ACTUAL: " + str(current_value))
    else:
      print("Test FAILED. Check theThe method to subtract two formulas with an exception because the second one contains more elements of a specific type than the first one. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))
  except IncorrectFormula as incfor:
    print("Test PASSED. The method to subtract two formulas where the second one contains a higher number of elements than the first one is working.")

  print("=================================================================.")
  print("Test Case 14: Testing the verification of a fragment mz from a formula with a null adduct")
  print("=================================================================.")
  try:
    formula_1 = 'C9H9O2S2'
    adduct = None 
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    current_value = my_formula_1.check_possible_fragment_mz(fragment_mz=200,ppm=50)
    expected_value = True
    if(current_value == expected_value):
      print("Test PASSED. The method to subtract two formulas where the second one contains a higher number of elements than the first one is working.")
    else:
      print("Test FAILED. Check theThe method to subtract two formulas with an exception because the second one contains more elements of a specific type than the first one. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check theThe method to subtract two formulas with an exception because the second one contains more elements of a specific type than the first one. Check the exception. ")


  print("=================================================================.")
  print("Test Case 15: Testing the verification of a fragment mz from a formula with an adduct that ")
  print("=================================================================.")
  try:
    formula_1 = 'C9H9O2S2'
    adduct = '[M-C2+H]+' 
    my_formula_1 = Formula.formula_from_str(formula_1, adduct)
    current_value = my_formula_1.check_possible_fragment_mz(fragment_mz=200,ppm=50)
    expected_value = False
    if(current_value == expected_value):
      print("Test PASSED. The method to subtract two formulas where the second one contains a higher number of elements than the first one is working.")
    else:
      print("Test FAILED. Check theThe method to subtract two formulas with an exception because the second one contains more elements of a specific type than the first one. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))
  except IncorrectFormula as incfor:
    print(incfor)
    print("Test FAILED. Check theThe method to subtract two formulas with an exception because the second one contains more elements of a specific type than the first one. Check the exception: ")

  
if __name__ == "__main__":
  main()