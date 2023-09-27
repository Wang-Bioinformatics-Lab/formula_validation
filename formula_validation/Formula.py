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

@version :  0.0.2, 27 September 2023
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
# libraries to make requests to chemcalc
import urllib3
import json
import numpy as np


from formula_validation.Element import Element_type, element_weights
from formula_validation.IncorrectFormula import IncorrectFormula
from formula_validation.NotFoundElement import NotFoundElement

class Formula:
  __electron_weight=0.00054858

  __connection_to_chemcalc_ws = urllib3.PoolManager()
  __ccurl = 'https://www.chemcalc.org/chemcalc/mf' 
  
  __default_ppm = 50


  """
  Methods:
    Constructor
      formula_from_str_hill(formula_str: str, adduct: str) -> 'Formula'
      formula_from_str(formula_str: str, adduct: str) -> 'Formula'
      formula_from_smiles(smiles: str, adduct: str) -> 'Formula'
      formula_from_inchi(inchi: str, adduct: str) -> 'Formula'
      __init__(self, elements: Dict[Union['Element_type',str], int], adduct: Union['Adduct', str], charge: int = 0, charge_type: str='', metadata: Dict=None): Initializes a Formula object with a dictionary of chemical elements and an optional adduct.
    Comparison and Representation
      __eq__(self, other): Checks if two Formula objects are equal.
      __str__(self): Returns a string representation of the Formula object.
      get_final_formula_with_adduct(self) -> str: Returns the final formula with the adduct as a string.
      __repr__(self): Returns a string representation of the Formula object.
      __hash__(self): Returns the hash value of the Formula object.
    Mathematical Operations
      __add__(self, other: 'Formula'): Adds another formula to the current one and keeps the adduct.
      __sub__(self, other: 'Formula'): Subtracts another formula from the current one and keeps the adduct.
      __mul__(self, num_to_multiply: int): Multiplies the formula by a specified number.
    Access and Information
      get_elements(self) -> Dict['Element_type',int]: Gets a copy of the dictionary of chemical elements and their counts in the formula.
      get_monoisotopic_mass(self) -> float: Gets the monoisotopic mass of the formula.
      get_monoisotopic_mass_with_adduct(self) -> float: Gets the monoisotopic mass of the formula, taking into account the adduct.
    Mass Comparison
      check_monoisotopic_mass(self, external_mass: Union[float,int], mass_tolerance_in_ppm: Union[int, float]) -> bool: Checks if the monoisotopic mass of the formula is within a specified mass tolerance of an external mass.
      check_monoisotopic_mass_with_adduct(self, external_mass: Union[float,int], mass_tolerance_in_ppm: Union[int, float]) -> bool: Checks if the monoisotopic mass of the formula, considering the adduct, is within a specified mass tolerance of an external mass.
      ppm_difference_with_exp_mass(self, reference_monoisotopic_mass: Union[float,int]) -> float: Calculates the ppm difference between the monoisotopic mass of the formula and a reference mass.
      absolute_to_ppm(reference_monoisotopic_mass: Union[float,int], mass_to_compare: Union[float,int]) -> float: Converts an absolute mass difference to ppm.
      ppm_to_absolute(reference_monoisotopic_mass: Union[float,int], ppm: Union[float,int]) -> float: Converts ppm to an absolute mass difference.
    Fragment Analysis
      check_possible_fragment_mz(self, fragment_mz: Union[float, int], ppm: Union[float, int]): Checks if a fragment mass can be explained by the formula and adduct.
      percentage_intensity_fragments_explained_by_formula(self, fragments_mz_intensities: Dict[Union[float, int], Union[float, int]], ppm: Union[float, int]): Calculates the percentage of intensity of fragments explained by the formula and adduct.
  Private Methods:
    __calculate_monoisotopic_mass(self) -> float: Calculates the monoisotopic mass of the formula.
    __calculate_monoisotopic_mass_with_adduct(self) -> float: Calculates the monoisotopic mass of the formula, taking into account the adduct.
  Static Methods: 
    Defines static methods to create Formula objects from different notations (Hill, SMILES, InChI).
"""
  
  def __init__(self, elements: Dict[Union['Element_type',str], int], adduct: Union['Adduct', str], charge: int = 0, charge_type: str='', metadata: Dict=None):
    """
    Constructor for the Formula class.
    Args:
      element_type (dict of Element_type and int): dictionary containing elements and their apps in a formula. If an element appears more than once, its appearences will be updated. Example {'C': 48, 'H': 72, 'N': 10, 'O': 12}
      adduct: string representing an adduct in the form '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+' where the charge is specified at the end . It should start with a [, then contain the multimer number followed by an M, then the adduct formula with a +-, the closing ], and the number of charges indicated by a number and the symbol +-
    
    Raises:
      IncorrectFormula: if the number of appearances is <=0
      NotFoundElement: if the dict contains elements that not a chemical element
      IncorrectAdduct: if the adduct is not a valid adduct with the format: '[M+C2H2O-H]-'
    
    Returns:
      Formula: a new instance of a formula
    """
    
    # The connection pool can be full so we to open a new one if not initialized
    if Formula.__connection_to_chemcalc_ws is None:
      try:
        Formula.__connection_to_chemcalc_ws = urllib3.PoolManager()
      except Exception as error:
        print("Error: Connection not established {}".format(error))
    
    self.metadata = metadata
    
    from collections.abc import Iterable
    from formula_validation.Adduct import Adduct

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
    self.__charge = charge
    
    if charge_type in ('','+','-'):
      self.__charge_type = charge_type
      self.__monoisotopic_mass = self.__calculate_monoisotopic_mass()
    else:
      raise IncorrectFormula("charge_type " + charge_type +  " invalid. It should be +, - or empty")
    
    
    if adduct is None:
      self.__adduct=None
    elif isinstance(adduct, str):
      if adduct == 'None' or adduct == '':
        self.__adduct=None
      else: 
        self.__adduct = Adduct(adduct)
    elif isinstance(adduct, Adduct): 
        self.__adduct = adduct
    else: 
      raise IncorrectFormula("The adduct " + str(adduct) + " is not a valid adduct with the format: '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+' ")
    
    self.__monoisotopic_mass_with_adduct = self.__calculate_monoisotopic_mass_with_adduct()

    
  def __eq__(self, other):
    """
      Check if two Formula objects are equal.

      Args:
        other (Formula): The other Formula object to compare with.

      Returns:
        bool: True if the two Formula objects are equal, otherwise False.
    """
    if isinstance(other, Formula):
      return self.__elements == other.get_elements() and self.__adduct == other.__adduct
    else:
      return False
    

  def __str__(self):
    """
      Return a string representation of the Formula object.

      Returns:
        str: A string representation of the Formula object in the format 'C4H5N6Na+[M+H]+'
    """
    formula_string = "".join( str(key.name) + (str(value) if value > 1 else "") for key,value in self.__elements.items())
    if self.__charge_type != '':
      charge_str = '' if self.__charge == 1 else str(self.__charge)
      formula_string = formula_string + self.__charge_type + charge_str
    adduct_str = "" if self.__adduct == None else "+" + str(self.__adduct)

    formula_string = formula_string + adduct_str
    return formula_string
  
  def get_final_formula_with_adduct(self) -> float:
    """
      Return a string representation of the final formula plus or minus de the adduct. 

      Returns:
        str: A string representation of the Formula object (C12H3N3O+[M-H2O+H]+) in the format '[C12H2N3]+'
    """
    if self.__adduct is None:
      return str(self)
    final_formula = self * self.__adduct.get_multimer()
    if self.__adduct.get_adduct_charge() == 0:
      formula_plus = self.__adduct.get_formula_plus()
      formula_minus = self.__adduct.get_formula_minus()
      final_formula = final_formula+formula_plus
      final_formula = final_formula-formula_minus
      return str(final_formula)
    else:
      if self.__charge_type == '+':
        own_charge = self.__charge 
      elif self.__charge_type == '-':
        own_charge = -self.__charge 
      else:
        own_charge = 0
      if self.__adduct.get_adduct_charge_type() == '+':
        adduct_charge = self.__adduct.get_adduct_charge() 
      elif self.__adduct.get_adduct_charge_type() == '-':
        adduct_charge = -self.__adduct.get_adduct_charge() 
      final_charge = own_charge + adduct_charge
      
      if final_charge == 0:
        return "".join( str(key.name) + (str(value) if value > 1 else "") for key,value in final_formula.get_elements().items())
      elif final_charge == 1:
        final_charge_str = "+" 
      elif final_charge > 1:
        final_charge_str = "+" + str(final_charge) 
      elif final_charge == -1:
        final_charge_str = "-" 
      else:
        final_charge_str = "-" + str(final_charge) 
            
      formula_plus = self.__adduct.get_formula_plus()
      formula_minus = self.__adduct.get_formula_minus()
      final_formula = final_formula+formula_plus
      final_formula = final_formula-formula_minus
      
      formula_string = "".join( str(key.name) + (str(value) if value > 1 else "") for key,value in final_formula.get_elements().items())
      formula_string = "[" + formula_string + "]" + final_charge_str
      return formula_string
    

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
      charge = -self.__charge if self.__charge_type == '-' else self.__charge
      charge = charge-other.__charge if other.__charge_type == '-' else charge+other.__charge
      if charge == 0:
        charge_type =''
      elif charge > 0:
        charge_type ='+'
      else: 
        charge_type ='-'
      return Formula(new_formula_dict, None, charge, charge_type)
    else:
        raise IncorrectFormula("other should be a formula and is a " + str(type(other)))
    

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
      # elements to remove are saved when the result is 0, since they should be later deleted. They cannot be deleted inmediately because the items are being iterated.
      elements_to_remove = set()
      for element, appearances in new_formula_dict.items():
        if appearances == 0:
          elements_to_remove.add(element)
        elif appearances < 0:
          raise IncorrectFormula("The subtraction of these two formulas contain a negative number of {element}")  
      for element_to_remove in elements_to_remove:
        del new_formula_dict[element_to_remove]
      
      charge = -self.__charge if self.__charge_type == '-' else self.__charge
      charge = charge+other.__charge if other.__charge_type == '-' else charge-other.__charge
      if charge == 0:
        charge_type =''
      elif charge > 0:
        charge_type ='+'
      else: 
        charge_type ='-'
      return Formula(new_formula_dict, None, charge, charge_type)
    
      
    else:
        raise IncorrectFormula("other should be a formula and is a " + str(type(other)))
    
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
      return Formula(new_formula_dict, self.__adduct, self.__charge, self.__charge_type)
    else:
        raise IncorrectFormula("other should be a formula and is a " + str(type(Formula)))
    
  @staticmethod
  def formula_from_str_hill(formula_str: str, adduct: str, metadata: dict=None) -> 'Formula':
    """
      Static method to create a Formula object from a chemical formula string in Hill notation.

      Args:
        formula_str (str): A string representing a molecular formula in Hill notation. Example: 'C4H5N6Na'. Other example 'C4H5N6Na+'
        adduct (str): A string representing an adduct in the form '[M+C2H2O-H]-', '[M-3H2O+2H]2+' or '[5M+Ca]2+' where the charge is specified at the end.
        metadata (dict): Optional argument to include a dict of metadata, defaults to None.

      Returns:
        Formula: A new instance of the Formula class with the elements specified in the string.

      Raises:
        IncorrectFormula: If the number of appearances is <=0 or if the formula contains elements that are not valid chemical elements.
    """
    import re
    
    if re.search(r'(?<![A-Z])[a-z]', formula_str):
      raise IncorrectFormula("The formula contains elements that are not chemical Elements")
    # check for any character that is not a capital letter, a lowercase letter, or a number
    if re.search(r'[^a-zA-Z0-9+-]', formula_str):
      raise IncorrectFormula("The formula contains parenthesis or brackets")

    #pattern = r'([A-Z][a-z]*)(\d*)([+-]?)(\d*)?'
    pattern = r'([A-Z][a-z]*)(\d*)'
    elements = {}
    
    for element, appearances in re.findall(pattern, formula_str):
      appearances = int(appearances) if appearances else 1
      elements[element] = elements.get(element, 0) + appearances
    
    charge_pattern = r'([-+]\d*)$'
    charge_match = re.search(charge_pattern, formula_str)

    if charge_match:
      charge_type = charge_match.group(1)[0]  # Capture the '+' or '-' symbol
      charge_value = charge_match.group(1)[1:]  # Capture the numeric part of the charge
      # If charge_value is empty, set charge to 1; otherwise, convert it to an integer
      charge = 1 if not charge_value else int(charge_value)
    else:
      charge = 0
      charge_type = ''

    return Formula(elements, adduct, charge, charge_type, metadata=metadata)
  
  @staticmethod
  def formula_from_str(formula_str: str, adduct: str, no_api: bool=False, metadata: bool=None) -> 'Formula':
    """
      Static method to create a Formula object from a chemical formula string.

      Args:
        formula_str (str): A string representing a molecular formula. Example: 'C4H5N6Na'.
        adduct (str): A string representing an adduct in the form '[M+C2H2O-H]-', '[M-3H2O+2H]2+' or '[5M+Ca]2+' where the charge is specified at the end.
        no_api (bool): Disables api calls for formula resolution.
        metadata (dict): Optional argument to include a dict of metadata, defaults to None.

      Returns:
        Formula: A new instance of the Formula class with the elements specified in the string.

      Raises:
        IncorrectFormula: If the number of appearances is <=0 or if the formula contains elements that are not valid chemical elements.
    """
    # First, we assume that formula is simple and it is not necessary to process it from the web server. 
    try:
      simple_formula = Formula.formula_from_str_hill(formula_str, adduct, metadata)
      return simple_formula
    except Exception as e:
      # If the formula could not be processed directly, then it is processed by chemcalc
      pass

    if not no_api:
      params = {'mf': formula_str,
        'isotopomers': 'jcamp,xy'
      }
      response = Formula.__connection_to_chemcalc_ws.request("GET",Formula.__ccurl, fields=params, retries = 3)

      # Read the output and convert it from JSON into a Python dictionary
      if(response.status != 200):
        raise IncorrectFormula("The formula " + formula_str + " was not parseable to a correct formula")
      
      data = json.loads(response.data)
      
      mf_hill = data['mf']
      return Formula.formula_from_str_hill(mf_hill, adduct, metadata)
    else:
      return None
    
  @staticmethod
  def formula_from_smiles(smiles: str, adduct: str, no_api: bool=False, metadata: Dict=None) -> 'Formula':
    """
      Static method to create a Formula object from a SMILES (Simplified Molecular Input Line Entry System) string.

      Args:
        smiles (str): A string representing a molecular structure in SMILES notation. Example: CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC
        adduct (str): A string representing an adduct in the form '[M+C2H2O-H]-', '[M-3H2O+2H]2+' or '[5M+Ca]2+' where the charge is specified at the end.
        no_api (bool): Disables api calls for formula resolution.
        metadata (dict): Optional argument to include a dict of metadata, defaults to None.

      Returns:
        Formula: A new instance of the Formula class according to the molecular structure.

      Raises:
        IncorrectFormula: If the SMILES string does not represent a valid molecular structure.
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
    return Formula.formula_from_str(formula, adduct, no_api, metadata=metadata)
    
  @staticmethod
  def formula_from_inchi(inchi: str, adduct: str, no_api:bool=False, metadata: Dict=None) -> 'Formula':
    """
      Static method to create a Formula object from an InChI (International Chemical Identifier) string.

      Args:
        inchi (str): A string representing a molecular structure in InChI notation. Example: InChI=1S/C45H73N5O10S3/c1-14-17-24(6)34(52)26(8)37-25(7)30(58-13)18-31-46-29(19-61-31)39-49-45(12,21-62-39)43-50-44(11,20-63-43)42(57)48-32(22(4)15-2)35(53)27(9)40(55)59-36(23(5)16-3)38(54)47-33(28(10)51)41(56)60-37/h19,22-28,30,32-37,51-53H,14-18,20-21H2,1-13H3,(H,47,54)(H,48,57)/t22-,23-,24+,25-,26-,27+,28+,30-,32-,33-,34-,35-,36-,37-,44+,45+/m0/s1
        adduct (str): A string representing an adduct in the form '[M+C2H2O-H]-', '[M-3H2O+2H]2+' or '[5M+Ca]2+' where the charge is specified at the end.
        no_api (bool): Disables api calls for formula resolution.
        metadata (dict): Optional argument to include a dict of metadata, defaults to None.

      Returns:
        Formula: A new instance of the Formula class according to the molecular structure.

      Raises:
        IncorrectFormula: If the InChI string does not represent a valid molecular structure.
    """

    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    if not inchi.startswith('InChI='):
      raise IncorrectFormula(inchi)
  
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
      raise IncorrectFormula(inchi)
    formula = CalcMolFormula(mol)
    return Formula.formula_from_str(formula, adduct, no_api, metadata)
      
  def get_elements(self) -> Dict['Element_type',int]:
    """
      Get a copy of the dictionary of chemical elements and their counts in the formula.

      Returns:
        Dict['Element_type', int]: A dictionary containing chemical elements as keys and their respective counts as values.
        """
    return self.__elements.copy()

  def __calculate_monoisotopic_mass(self) -> float:

    monoisotopic_mass = 0
    for element, appearances in self.__elements.items():
      monoisotopic_mass += element_weights[element] * appearances

    if self.__charge_type=='+':
      electrons_weight = -Formula.__electron_weight*self.__charge
    elif self.__charge_type=='-':
      electrons_weight = Formula.__electron_weight*self.__charge
    elif self.__charge_type=='':
      electrons_weight = 0
    else:
      raise IncorrectFormula("The formula contains a wrong charge type")
    
    monoisotopic_mass += electrons_weight
    adduct_charge_to_divide = self.__charge if self.__charge != 0 else 1
    monoisotopic_mass = monoisotopic_mass / adduct_charge_to_divide

    return monoisotopic_mass
        
  def get_monoisotopic_mass(self) -> float:
    """
        Get the monoisotopic mass of the formula, taking into account the adduct.

        Returns:
          float: The monoisotopic mass of the formula with the adduct considered.
    """

    return self.__monoisotopic_mass
  
  def __calculate_monoisotopic_mass_with_adduct(self) -> float:
    """
      Get the monoisotopic mass of the formula, taking into account the adduct.

      Returns:
        float: The monoisotopic mass of the formula with the adduct considered.
    """
    monoisotopic_mass_with_adduct= self.get_monoisotopic_mass()
    if self.__adduct == None:
      return monoisotopic_mass_with_adduct
    
    partial_elements = {element: count * self.__adduct.get_multimer() for element, count in self.__elements.items()}
    if self.__charge_type=='+':
        partial_charge = self.__charge
    else:
      partial_charge = -self.__charge 
    
    if self.__adduct.get_adduct_charge_type()=='+':
      final_charge=partial_charge+self.__adduct.get_adduct_charge()
    elif self.__adduct.get_adduct_charge_type()=='-':
      final_charge=partial_charge-self.__adduct.get_adduct_charge()
    elif self.__adduct.get_adduct_charge_type()=='':
      final_charge = partial_charge
    else:
      raise IncorrectFormula("The formula contains a wrong adduct")
    
    # Calculate the final number of elements with the multimers and the adduct
    monoisotopic_mass_with_adduct = 0
    formula_plus = self.__adduct.get_formula_plus()
    formula_minus = self.__adduct.get_formula_minus()
    
    for element, appearances in formula_plus.get_elements().items():
      if element in partial_elements:
        partial_elements[element] +=  appearances 
      else:
        partial_elements[element] =  appearances 
    
    for element, appearances in formula_minus.get_elements().items():
      if element in partial_elements:
        partial_elements[element] -= appearances 
      else:
        partial_elements[element] = appearances 
      if(partial_elements[element] < 0):
        raise IncorrectFormula("The formula contains a wrong adduct because the element " + str(element) + " is negative " + str(partial_elements[element]))

    for element, appearances in partial_elements.items():
      monoisotopic_mass_with_adduct += element_weights[element] * appearances
    

    if final_charge != 0:
      electrons_weight = -Formula.__electron_weight*final_charge
    else:
      electrons_weight = 0
    
    monoisotopic_mass_with_adduct += electrons_weight
    adduct_charge_to_divide = final_charge if final_charge != 0 else 1
    monoisotopic_mass_with_adduct = monoisotopic_mass_with_adduct / abs(adduct_charge_to_divide)
    
    return monoisotopic_mass_with_adduct

  def get_monoisotopic_mass_with_adduct(self) -> float:
    """
      Get the monoisotopic mass of the formula, taking into account the adduct.

      Returns:
        float: The monoisotopic mass of the formula with the adduct considered.
    """
    return self.__monoisotopic_mass_with_adduct
    

  def check_monoisotopic_mass(self, external_mass: Union[float,int], mass_tolerance_in_ppm: Union[int, float] = __default_ppm) -> bool:
    """
      Check if the monoisotopic mass of the formula is within a specified mass tolerance of an external mass.

      Args:
        external_mass (Union[float,int]): The external monoisotopic mass to compare with the formula's mass.
        mass_tolerance_in_ppm (Union[int, float]): The mass tolerance in parts per million (ppm) for the comparison.

      Returns:
        bool: True if the external mass is within the specified tolerance of the formula's mass, otherwise False.
    """
    import math
    abs_value_delta = Formula.ppm_to_absolute(self.get_monoisotopic_mass(), mass_tolerance_in_ppm)
    if math.isclose(self.get_monoisotopic_mass(),external_mass, abs_tol=abs_value_delta):
      return True
    else:
      return False
  
  def check_monoisotopic_mass_with_adduct(self, external_mass: Union[float,int], mass_tolerance_in_ppm: Union[int, float] = __default_ppm) -> bool:
    """
      Check if the monoisotopic mass of the formula, considering the adduct, is within a specified mass tolerance of an external mass.

      Args:
        external_mass (Union[float,int]): The external monoisotopic mass to compare with the formula's mass with the adduct.
        mass_tolerance_in_ppm (Union[int, float]): The mass tolerance in parts per million (ppm) for the comparison.

      Returns:
        bool: True if the external mass is within the specified tolerance of the formula's mass with the adduct, otherwise False.
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
      
      Raise: 
        A exception if reference_monoisotopic_mass or ppm are not numbers  
      
      Returns: 
        The ppms between the monoisotopic mass of the formula taking into account the adduct and the experimental mass detected
      
    """
    return Formula.absolute_to_ppm(self.get_monoisotopic_mass_with_adduct(), reference_monoisotopic_mass)

  @staticmethod
  def absolute_to_ppm(reference_monoisotopic_mass: Union[float,int], mass_to_compare: Union[float,int]) -> float:
    """
      Args:
        reference_monoisotopic_mass (numeric): monoisotopic mass of reference
        mass_to_compare(numeric): mass to compare
      
      Raise: 
        A exception if reference_monoisotopic_mass or ppm are not numbers
        
      Returns: 
        The ppms between the reference_monoisotopic_mass mass and mass_to_compare
    
    """
    
    ppm_diff = abs((reference_monoisotopic_mass - mass_to_compare) / reference_monoisotopic_mass) * 1000000.0
    return ppm_diff

  @staticmethod
  def ppm_to_absolute(reference_monoisotopic_mass: Union[float,int], ppm: Union[float,int] = __default_ppm) -> float:
    """
      Args:
        reference_monoisotopic_mass (numeric): monoisotopic mass of reference
        ppm (numeric): ppm of the reference monoisotopic mass

      Raise: 
        a exception if reference_monoisotopic_mass or ppm are not numbers

      Returns: 
        the absolute value of the ppm calculated
    """
    return (reference_monoisotopic_mass / 1000000.0) * ppm
  
  def check_possible_fragment_mz(self, fragment_mz: Union[float, int], ppm: Union[float, int] = __default_ppm):
    """
      Args:
        self (Formula): formula to check if the fragment could be generated by this formula according to the substructure and the adduct. If the adduct has a number of charges > 1, it will check the adducts for double and single charge
        fragment_mz (numeric): mz_mass to check if it can come from the structure
        ppm (numeric): ppm of the reference monoisotopic mass to check the potential formulas
      
      Raise: 
        A exception if fragment_mz or ppm are not numbers
      
      Returns: 
        A boolean specifying if the fragment m_z can be explained from the formula and the adduct
    """

    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.rinterface_lib.sexp import NULLType

    robjects.r('''library(MassTools)''')

    mz = 200
    charge = 1
    
    top = 5
    final_formula = self
    final_charge = - self.__charge if self.__charge_type=='-' else self.__charge
    if self.__adduct is not None:
      final_charge= final_charge + self.__adduct.get_adduct_charge() if self.__adduct.get_adduct_charge_type() =='+' else final_charge - self.__adduct.get_adduct_charge()
      final_formula = final_formula + self.__adduct.get_formula_plus() if self.__adduct.get_formula_plus() is not None else final_formula
      final_formula = final_formula - self.__adduct.get_formula_minus() if self.__adduct.get_formula_minus() is not None else final_formula
    elements = {key: value for key, value in final_formula.get_elements().items() if value > 0}
    

    min_elements = "".join([f"{element.name}0" for element, count in elements.items()])
    max_elements = "".join([f"{element.name}{count}" for element, count in elements.items()])
    
    elements_list=[key.name for key in elements]
    # Convert the filter strings to R format
    robjects.r(f'minElements = "{min_elements}"')
    robjects.r(f'maxElements = "{max_elements}"')

    r_vector_str = '", "'.join(elements_list)

    # Create the R code to initialize the 'elements' list
    elements = robjects.r(f'elements = Rdisop::initializeElements(c("{r_vector_str}"))')


    #robjects.r('elements = Rdisop::initializeElements(c("C", "H", "N", "O", "S"))')



    # Create the filter list in R format
    filter_r_code = '''
    list(
        minElements = minElements,
        maxElements = maxElements,
        maxCounts = TRUE,
        HCratio = FALSE,
        moreRatios = TRUE
    )
    '''

    # Run the R code in Python
    robjects.r('mz <- {}'.format(fragment_mz))
    robjects.r('z <- {}'.format(final_charge))
    robjects.r('ppm <- {}'.format(ppm))
    robjects.r('top <- {}'.format(top))
    robjects.r('mfs <- calcMF(mz = mz, z = z, ppm = ppm, top = top, elements = elements, Filters = {})'.format(filter_r_code))
    # This is actually an R dataframe with information about the formula that explains the peak. TO DO: use in future
    mfs = robjects.r('mfs')
    if isinstance(mfs, NULLType):
      return False
    return True
  
  def percentage_intensity_fragments_explained_by_formula(self, fragments_mz_intensities: Dict[Union[float, int], Union[float, int]], ppm: Union[float, int] = __default_ppm):
    """
      Args:
        self (Formula): formula to check if the fragment could be generated by this formula according to the substructure and the adduct. If the adduct has a number of charges > 1, it will check the adducts for double and single charge
        fragments_mz_intensities (dict[numeric,numeric]): dict containing fragments in the format m/z, intensity. 
        ppm (numeric): ppm of the reference monoisotopic mass to check the potential formulas
      
      Raise: 
        a exception if fragment_mz or ppm are not numbers
      
      Returns: 
        the percentage of intensity of fragments explained according to the formula and adduct
      
    """
    
    total_fragment_intensities = sum(value for value in fragments_mz_intensities.values() if isinstance(value, (int, float)))
    if total_fragment_intensities <= 1e-16:
      # We can't guarentee any names so we'll add all possbily useful information
      print(f"Zero intensity found for entry with formula {self.__elements} and metadata {self.metadata}")
      return np.nan
    
    explained_intensities = 0
    for mz, intensity in fragments_mz_intensities.items():
      if self.check_possible_fragment_mz(mz,ppm):
        explained_intensities +=intensity
    return explained_intensities / total_fragment_intensities


if __name__ == '__main__':
  smiles_1 = 'CCCCCCC[C@@H](C/C=C/CCC(=O)NC/C(=C/Cl)/[C@@]12[C@@H](O1)[C@H](CCC2=O)O)OC'
  adduct = '[M+C2H2O-H]-' 
  my_formula = Formula.formula_from_smiles(smiles_1, adduct)
  expected_value =496.24713858
  current_value = my_formula.get_monoisotopic_mass_with_adduct()