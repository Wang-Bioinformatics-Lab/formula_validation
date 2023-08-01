#!/usr/bin/env python

"""
This Python module contains not only the class Formula, but also the test of
this Formula class.

@contents :  This Python module contains not only the class Formula, but also the test of
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

from typing import Union, Dict
# 
import urllib3
import json

from Element import Element_type, element_weights
from IncorrectFormula import IncorrectFormula
from NotFoundElement import NotFoundElement


class Formula:
  __electron_weight=0.00054858

  __connection_to_chemcalc_ws = urllib3.PoolManager()
  __ccurl = 'https://www.chemcalc.org/chemcalc/mf' 
  


  """
  Methods:

  constructor(elements,adduct): (receives a dict of chemical elements and its apparenaces > 0) and a string representing an adduct in the form '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+' where the charge is specified at the end . It should start with a [, then contain the multimer number followed by an M, then the adduct formula with a +-, the closing ], and the number of charges indicated by a number and the symbol +-
  get_formula_from_str(formula_str): STATIC. Returns a formula from a string representing a formula
  get_formula_from_smiles(smiles): STATIC. Returns a formula from a SMILES representing a structure
  check_monoisotopic_mass(external_mass, mass_tolerance_in_ppm=50): check if the mass of the formula and a external mass have a difference higher than the mass_tolerance in ppm established
  get_monoisotopic_mass(): returns the mass of the formula
  get_monoisotopic_mass_with_adduct(): returns the monoisotopic mass of the formula +- the adduct
    """
  
  def __init__(self, elements: Dict[Union['Element_type',str], int], adduct: Union['Adduct', str]):
    """
    Args:
      element_type (dict of Element_type and int): dictionary containing elements and their apps in a formula. If an element appears more than once, its appearences will be updated. Example {'C': 48, 'H': 72, 'N': 10, 'O': 12}
      adduct: string representing an adduct in the form '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+' where the charge is specified at the end . It should start with a [, then contain the multimer number followed by an M, then the adduct formula with a +-, the closing ], and the number of charges indicated by a number and the symbol +-
    Returns:
      Formula: a new instance of a formula
    Raises:
      IncorrectFormula: if the number of appearances is <=0
      NotFoundElement: if the dict contains elements that not a chemical element
      IncorrectAdduct: if the adduct is not a valid adduct with the format: '[M+C2H2O-H]-'
    """
    
    # The connection pool can be full so we to open a new one if not initialized
    if Formula.__connection_to_chemcalc_ws is None:
      try:
        Formula.__connection_to_chemcalc_ws = urllib3.PoolManager()
      except Exception as error:
        print("Error: Connection not established {}".format(error))
    
    from collections.abc import Iterable
    from Adduct import Adduct
    self.__elements={}
    if isinstance(elements, dict):
      for element, appearances in elements.items():
        if not isinstance(appearances,int):
          raise IncorrectFormula(elements)
        elif appearances <= 0:
          raise IncorrectFormula(elements)
        if isinstance(element,Element_type):
          self.__elements[element] = self.__elements.get(element, 0) + appearances
        elif isinstance(element,str):
          try:
            element = Element_type[element]
            self.__elements[element] = self.__elements.get(element, 0) + appearances
          except KeyError as ke:
            raise NotFoundElement(element)
        else:
          raise IncorrectFormula(elements)
    else:
      raise IncorrectFormula(elements)
    self.__monoisotopic_mass = self.__calculate_monoisotopic_mass()
    
    if adduct is None:
      self.__adduct=None
    elif isinstance(adduct, str):
      if adduct == 'None' or adduct == '':
        self.__adduct=None
      else: 
        self.__adduct = Adduct(adduct)
    elif isinstance(adduct, Adduct): 
        self.__adduct = Adduct(adduct)
    else: 
      raise IncorrectFormula("The adduct " + str(adduct) + " is not a valid adduct with the format: '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+' ")
    self.__monoisotopic_mass_with_adduct = self.__calculate_monoisotopic_mass_with_adduct()

    
  def __eq__(self, other):
    if not isinstance(other, Formula):
      return False
    return self.__elements == other.get_elements()

  def __str__(self):
    
    return str(self.__elements)
  
  def __repr__(self):
    return str(self)

  def __hash__(self):
    return hash(frozenset(self.__elements.items()))
  
  def __add__(self, other: 'Formula'):
    """
    Args:
        other (Formula): another formula to add the elements with the current one and keeps the adduct of the current formula

    Raises:
        IncorrectFormula: if the object is not a formula

    Requirements: 
        module copy is used to make a deep copy of an object
    Returns:
        formula (Formula): a new formula that is the addition of the chemical elements from both formulas
    """
    import copy
    if isinstance(other, Formula):
      new_formula_dict = copy.deepcopy(self.__elements)
      for element,counts_in_other in other.__elements.items():
        new_formula_dict[element] = new_formula_dict.get(element,0) + counts_in_other
      return Formula(new_formula_dict, self.__adduct)
    else:
        raise IncorrectFormula("other should be a formula and is a " + type(Formula))
    

  def __sub__(self, other: 'Formula'):
    """
    Args:
        other (Formula): another formula to subtract the elements from the current one and keeps the adduct of the current formula

    Raises:
        IncorrectFormula: if the object is not a formula

    Requirements: 
        module copy is used to make a deep copy of an object
    Returns:
        formula (Formula): a new formula that is the subtraction of the chemical elements 
    """
    import copy
    if isinstance(other, Formula):
      new_formula_dict = copy.deepcopy(self.__elements)
      for element,counts_in_other in other.__elements.items():
        new_formula_dict[element] = new_formula_dict.get(element,0) - counts_in_other
      # I save the elements to remove when the result is 0
      elements_to_remove = set()
      for element, appearances in new_formula_dict.items():
        if appearances == 0:
          elements_to_remove.add(element)
        elif appearances < 0:
          raise IncorrectFormula("The addition of these two formulas contain a negative number of {element}")  
      for element_to_remove in elements_to_remove:
        del new_formula_dict[element_to_remove]
      return Formula(new_formula_dict, self.__adduct)
      
    else:
        raise IncorrectFormula("other should be a formula and is a " + type(Formula))
    
  def __mul__(self, num_to_multiply: int):
    """
    Args:
        other (Formula): another formula to add the elements with the current one

    Raises:
        IncorrectFormula: if the object is not a formula

    Requirements: 
        module copy is used to make a deep copy of an object
    Returns:
        formula (Formula): a new formula that is the addition of the chemical elements from both formulas
    """
    import copy
    if isinstance(num_to_multiply, int):

      new_formula_dict = copy.deepcopy(self.__elements)
      for element,counts_in_other in self.__elements.items():
        new_formula_dict[element] = new_formula_dict.get(element,0) * num_to_multiply
      return Formula(new_formula_dict, self.__adduct)
    else:
        raise IncorrectFormula("other should be a formula and is a " + type(Formula))
    

  def formula_from_str_hill(formula_str: str, adduct: str) -> 'Formula':
    """
    Args:
      formula_str (str): represents a molecular formula as a string of type [Element_type][NumberOfOccurences]: C4H5N6Na. It can contain parenthesis. 
      adduct (str): adduct representing the adduct formed by the molecular formula expressed by the form '[M+C2H2O-H]-'
    Returns:
      Formula: a new instance of a formula with the elements specified in the string
    Raises:
      IncorrectFormula: if the number of appearances is <=0
      NotFoundElement: if the dict contains elements that not a chemical element
    """
    import re
    pattern = r'([A-Z][a-z]*)(\d*)'
    elements = {}
    
    for element, appearances in re.findall(pattern, formula_str):
      appearances = int(appearances) if appearances else 1
      elements[element] = elements.get(element, 0) + appearances
      
    return Formula(elements, adduct)
  
  def formula_from_str(formula_str: str, adduct: str) -> 'Formula':
    """
    Args:
      formula_str (str): represents a molecular formula as a string of type [Element_type][NumberOfOccurences]: C4H5N6Na. It can contain parenthesis. 
      adduct (str): adduct representing the adduct formed by the molecular formula expressed by the form '[M+C2H2O-H]-'
    Returns:
      Formula: a new instance of a formula with the elements specified in the string
    Raises:
      IncorrectFormula: if the number of appearances is <=0
      NotFoundElement: if the dict contains elements that not a chemical element
    """
    # First, we assume that formula is simple and it is not necessary to process it from the web server. 
    try:
      simple_formula = Formula.formula_from_str_hill(mf_hill, adduct)
      return simple_formula
    except exception as e:
      # If the formula could not be processed directly, then it is processed by chemcalc
      pass

    params = {'mf': formula_str,
      'isotopomers': 'jcamp,xy'
    }
    response = Formula.__connection_to_chemcalc_ws.request("GET",Formula.__ccurl, fields=params, retries = 3)

    # Read the output and convert it from JSON into a Python dictionary
    if(response.status != 200):
      raise IncorrectFormula("The formula " + formula_str + " was not parseable to a correct formula")
    
    data = response.json()
    
    mf_hill = data['mf']
    return Formula.formula_from_str_hill(mf_hill, adduct)
    
    
  def formula_from_smiles(smiles: str, adduct: str) -> 'Formula':
    """
    Args:
      smiles (str): represents a molecular structure as a string. Example: CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC
      adduct (str): adduct representing the adduct formed by the molecular formula expressed by the form '[M+C2H2O-H]-'
    Returns:
      Formula: according to the structure
    Raises:
      IncorrectFormula: if the SMILES does not represent a structure
    """
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    if smiles =='' or smiles=='nan' or smiles == 'None' or smiles == None:
      raise IncorrectFormula(smiles)
    elif smiles.startswith('InChI='):
      return Formula.formula_from_inchi(smiles,adduct)
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
      raise IncorrectFormula(smiles)
    formula = CalcMolFormula(mol)
    return Formula.formula_from_str(formula, adduct)
    
    
  def formula_from_inchi(inchi: str, adduct: str) -> 'Formula':
    """
    Args:
      inchi (str): represents a molecular structure as a string. Example: InChI=1S/C45H73N5O10S3/c1-14-17-24(6)34(52)26(8)37-25(7)30(58-13)18-31-46-29(19-61-31)39-49-45(12,21-62-39)43-50-44(11,20-63-43)42(57)48-32(22(4)15-2)35(53)27(9)40(55)59-36(23(5)16-3)38(54)47-33(28(10)51)41(56)60-37/h19,22-28,30,32-37,51-53H,14-18,20-21H2,1-13H3,(H,47,54)(H,48,57)/t22-,23-,24+,25-,26-,27+,28+,30-,32-,33-,34-,35-,36-,37-,44+,45+/m0/s1
      adduct (str): adduct representing the adduct formed by the molecular formula expressed by the form '[M+C2H2O-H]-'
    Returns:
      Formula: according to the structure
    Raises:
      IncorrectFormula: if the SMILES does not represent a structure
    """
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    if not inchi.startswith('InChI='):
      raise IncorrectFormula(inchi)
  
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
      raise IncorrectFormula(inchi)
    formula = CalcMolFormula(mol)
    return Formula.formula_from_str(formula, adduct)
      
  def get_elements(self) -> Dict['Element_type',int]:
    """
    Returns: A copy of the dictionary of the elements so the formula cannot be mutated
    """
    return self.__elements.copy()

  def __calculate_monoisotopic_mass(self) -> float:

    monoisotopic_mass = 0
    for element, appearances in self.__elements.items():
      monoisotopic_mass += element_weights[element] * appearances
    return monoisotopic_mass
        
  def get_monoisotopic_mass(self) -> float:
    """
    Returns: the monoisotopic mass of the formula taking into account the adduct formed, such as '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+'
    """
    return self.__monoisotopic_mass
  
  def __calculate_monoisotopic_mass_with_adduct(self) -> float:
    """
    Returns: the monoisotopic mass of the formula taking into account the adduct formed, such as '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+'
    """
    monoisotopic_mass_with_adduct= self.get_monoisotopic_mass()
    if self.__adduct == None:
      return monoisotopic_mass_with_adduct
    if self.__adduct.get_multimer()>1:
      monoisotopic_mass_with_adduct = monoisotopic_mass_with_adduct * self.__adduct.get_multimer()

    if self.__adduct.get_adduct_charge_type()=='+':
      electrons_weight = -Formula.__electron_weight*self.__adduct.get_adduct_charge()
    elif self.__adduct.get_adduct_charge_type()=='-':
      electrons_weight = Formula.__electron_weight*self.__adduct.get_adduct_charge()
    elif self.__adduct.get_adduct_charge_type()=='':
      electrons_weight = 0
    else:
      raise IncorrectFormula("The formula contains a wrong adduct")
    
    monoisotopic_mass_with_adduct += electrons_weight
    monoisotopic_mass_with_adduct += self.__adduct.get_adduct_mass()

    monoisotopic_mass_with_adduct = monoisotopic_mass_with_adduct / self.__adduct.get_adduct_charge()
    
    return monoisotopic_mass_with_adduct

  def get_monoisotopic_mass_with_adduct(self) -> float:
    """
    Returns: the monoisotopic mass of the formula taking into account the adduct coupled to the structure
    """
    return self.__monoisotopic_mass_with_adduct
    

  def check_monoisotopic_mass(self, external_mass: Union[float,int], mass_tolerance_in_ppm: Union[int, float] =50) -> bool:
    """
    Args:
      external_mass (numeric): represents a monoisotopic mass to be compared with the mass of the formula
      mass_tolerance_in_ppm (numeric): mass tolerance permitted
    Returns: wether the external_mass is within the mass of the formula +- the tolerance established in ppm
    Raise: a exception if external_mass or mass_tolerance_in_ppm are not numbers
    """
    import math
    abs_value_delta = Formula.ppm_to_absolute(self.get_monoisotopic_mass(), mass_tolerance_in_ppm)
    if math.isclose(self.get_monoisotopic_mass(),external_mass, abs_tol=abs_value_delta):
      return True
    else:
      return False
  
  def check_monoisotopic_mass_with_adduct(self, external_mass: Union[float,int], mass_tolerance_in_ppm: Union[int, float] =50) -> bool:
    """
    Args:
      external_mass (numeric): represents a monoisotopic mass to be compared with the mass of the formula
      mass_tolerance_in_ppm (numeric): mass tolerance permitted
    Returns: wether the external_mass is within the mass of the formula +- the tolerance established in ppm
    Raise: a exception if external_mass or mass_tolerance_in_ppm are not numbers
    """
    import math
    abs_value_delta = Formula.ppm_to_absolute(self.get_monoisotopic_mass_with_adduct(), mass_tolerance_in_ppm)
    if math.isclose(self.get_monoisotopic_mass_with_adduct(),external_mass, abs_tol=abs_value_delta):
      return True
    else:
      return False
  

  def ppm_difference_with_exp_mass(self, reference_monoisotopic_mass: Union[float,int]) -> float:
    """
    Args:
      reference_monoisotopic_mass (numeric): monoisotopic mass of reference
      
    Returns: the ppms between the monoisotopic mass of the formula taking into account the adduct and the experimental mass detected
    Raise: a exception if reference_monoisotopic_mass or ppm are not numbers
    """
    print(self.get_monoisotopic_mass_with_adduct(), reference_monoisotopic_mass)
    return Formula.absolute_to_ppm(self.get_monoisotopic_mass_with_adduct(), reference_monoisotopic_mass)

  
  def absolute_to_ppm(reference_monoisotopic_mass: Union[float,int], mass_to_compare: Union[float,int]) -> float:
    """
    Args:
      reference_monoisotopic_mass (numeric): monoisotopic mass of reference
      mass_to_compare(numeric): mass to compare
    Returns: the ppms between the reference_monoisotopic_mass mass and mass_to_compare
    Raise: a exception if reference_monoisotopic_mass or ppm are not numbers
    """
    
    ppm_diff = abs((reference_monoisotopic_mass - mass_to_compare) / reference_monoisotopic_mass) * 1000000.0
    return ppm_diff


  def ppm_to_absolute(reference_monoisotopic_mass: Union[float,int], ppm: Union[float,int] = 50) -> float:
    """
    Args:
      reference_monoisotopic_mass (numeric): monoisotopic mass of reference
      ppm (numeric): ppm of the reference monoisotopic mass
    Returns: the absolute value of the ppm calculated
    Raise: a exception if reference_monoisotopic_mass or ppm are not numbers
    """
    return (reference_monoisotopic_mass / 1000000.0) * ppm
  

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
    print("Test FAILED. Check the constructor of formula from a string")

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
    print("Test FAILED. Check the constructor of formula from a string")

  
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
    print("Test FAILED. Check the constructor of formula from a string")

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
    print("Test FAILED. Check the function to calculate absolute values from ppms")

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
    print("Test FAILED. Check the function to check a formula and an external mass have a lower difference than the ppm establised. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))

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
    print("Test FAILED. Check the function to calculate the monoisotopic mass with an adduct in a single negative charge. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))
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
    print("Test FAILED. Check the function to calculate the monoisotopic mass with an adduct in a single positive charge. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))

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
    print("Test FAILED. Check the function to calculate the monoisotopic mass with an adduct in a single positive charge and multimer 5. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))

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
    print("Test FAILED. Check the function Calculating the ppm difference between a formula and an experimental value. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))

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
    print("Test FAILED. Check the method The function to add two formulas. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))

  print("=================================================================.")
  print("Test Case 12: Testing the subtraction of two formulas with a correct formula")
  print("=================================================================.")
  try:
    formula_1 = 'C4H4O2'
    formula_2 = 'C4H2O1'
    adduct = None # adduct weight = 18.01056
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
    print("Test FAILED. The method to subtract two formulas with proper result. Expected: " + str(expected_value) + " ACTUAL: " + str(current_value))


  print("=================================================================.")
  print("Test Case 13: Testing the subtraction of two formulas with netagive elements as a result")
  print("=================================================================.")
  try:
    formula_1 = 'C4H4O2'
    formula_2 = 'C4H6O1'
    adduct = None # adduct weight = 41.00328858
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



if __name__ == "__main__":
  main()