#!/usr/bin/env python

"""
This Python module contains not only the class Adduct, but also the test of
this Adduct class.

@contents :  This Python module contains not only the class Adduct, but also the test of
this Adduct class.
@project :  N/A
@program :  N/A
@file :  Adduct.py
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
__author__      = "Alberto Gil de la Fuente"
__copyright__   = "GPL License version 3"
from formula_validation.IncorrectAdduct import IncorrectAdduct
from formula_validation.Formula import Formula
import re

class Adduct:
  """
  Methods:
    get_multimer(self) -> int: Returns the number of multimers formed.
    get_formula_str(self) -> str: Returns the formula string representing the adduct (elements to add and subtract).
    get_formula_plus(self) -> 'Formula': Returns the Formula object representing the elements to add in the adduct.
    get_formula_minus(self) -> 'Formula': Returns the Formula object representing the elements to subtract in the adduct.
    get_adduct_mass(self) -> float: Returns the mass of the adduct.
    get_adduct_charge(self) -> int: Returns the number of charges (an integer greater than 0).
    get_adduct_charge_type(self) -> str: Returns the type of charge, either "+" or "-".
    Private Methods:

    __init__(self, adduct: str): Constructor for the Adduct class, parses the input adduct string.
    __eq__(self, other): Checks if two Adduct objects are equal.
    __str__(self): Returns a string representation of the Adduct object.
    __repr__(self): Returns a string representation of the Adduct object for debugging.
    __hash__(self): Returns a hash value for the Adduct object.
    __calculate_adduct_formula_to_add_and_subtract(self, adduct_formula: str) -> float: Calculates the formula to add and subtract elements based on the input adduct formula.
  """ 
  __electron_weight=0.00054858

  def __init__(self, adduct: str):
    
    """
      Args:
        adduct (str): String representing an adduct in the form '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+' where the charge is specified at the end . It should start with a [, then contain the multimer number followed by an M, then the adduct formula with a +-, the closing ], and the number of charges indicated by a number and the symbol +-
        It can contain as much subformulas as wanted.
      Properties:
        multimer (int): number of M in the structure where the adduct is present
        formula_plus (Formula): formula of the adduct to add to the mass. It will be None if they dont apply or if adduct is [xM], [xM]+ or [xM]- being x an integer > 1
        formula_minus (Formula): formula of the adduct to subtract from the mass. It will be None if they dont apply or if adduct is [xM], [xM]+ or [xM]- being x an integer > 1
        adduct_mass: mass to add or remove from the structure the adduct is present
        charges: number of charges
        charge_type (+ or -): if the adduct is negatively or positively charged
      Raises:
        IncorrectAdduct: if the adduct does not follow the format specified
    """
    
    pattern = r'\[(\d)*M([\+-].*?)\](\d)*([\+-])?'

    match = re.match(pattern, adduct)
    if match:
        self.__multimer = int(match.group(1).strip()) if match.group(1) is not None else 1
        self.__original_formula = match.group(2).strip()
        
        # self.__formula_plus, self.__formula_minus, self.__adduct_mass. One or both of them will be None if they dont apply or if adduct is [xM], [xM]+ or [xM]- being x an integer > 1
        self.__formula_plus, self.__formula_minus, self.__adduct_mass = self.__calculate_adduct_formula_to_add_and_subtract(self.__original_formula)
        
        if match.group(4) is not None:
          self.__charge = int(match.group(3).strip()) if match.group(3) is not None else 1
          self.__charge_type = match.group(4).strip() 
        else:
          self.__charge = 0
          self.__charge_type = ""

    else:
        raise IncorrectAdduct(adduct)
    
  def __eq__(self, other):
    if not isinstance(other, Adduct):
      return False
    return (self.__multimer == other.__multimer and self.__formula_plus == other.__formula_plus and self.__formula_minus == other.__formula_minus and self.__charge == other.__charge and self.__charge_type == other.__charge_type)

  def __str__(self):
    multimer = "" if self.__multimer == 1 else str(self.__multimer)
    charge = "" if self.__charge in(0,1) else str(self.__charge)
    
    return f"[{multimer}M{self.get_formula_str()}]{charge}{self.__charge_type}"

  
  def __repr__(self):
    return str(self)

  def __hash__(self):
    return hash(self.__multimer, self.__formula_plus, self.__formula_minus, self.__charge, self.__charge_type)


  def __calculate_adduct_formula_to_add_and_subtract(self, adduct_formula: str) -> float:
    """
      Args:
        formula (str): String representing the formula within an adduct in the form +HCOOH-H, +Ca, +H, +CH3COOH-H, etc..
      Returns:
        The formula with the elements of the adduct to add and other formula with the elements of the adduct to subtract. It unifies any type of formula with the total number of elements to add or subtract.
      Raises:
        IncorrectAdduct: if the adduct does not follow the format specified
    """

    # It will not use the dict from the formula, instead it creates a local one. This is done because the formula dict does not allow a negative number of elements, but the calculation of the formula of the adduct can have temporary and final states with negative number 
    # and it may have the need to create two different formulas to add some elements and subtract some other elements, like for example an adduct that is [M-2H+K]- that should subtract 2 hydrogens and add a potasium

    self.__adduct_mass=0
    pattern = r'([\+-])(\d)*([A-Za-z0-9]+)*'
    formulas = re.findall(pattern, adduct_formula)
    final_adduct_elements_dict = {}
    for symbol, number_subformulas, subformula_str in formulas:
      number_subformulas = int(number_subformulas) if number_subformulas else 1
      subformula=Formula.formula_from_str(subformula_str, 'None')
      dict_subformula=subformula.get_elements()
      
      if symbol == '+':
        for element,count in dict_subformula.items():
          number_elements = final_adduct_elements_dict[element] if element in final_adduct_elements_dict else 0
          final_adduct_elements_dict[element] = number_elements + count * number_subformulas
      elif symbol == '-': 
        for element,count in dict_subformula.items():
          number_elements = final_adduct_elements_dict[element] if element in final_adduct_elements_dict else 0
          final_adduct_elements_dict[element] = number_elements - count * number_subformulas
      else: 
        raise IncorrectAdduct(adduct_formula)
    
    elements_dict_to_add={}
    elements_dict_to_subtract={}
    for element, count in final_adduct_elements_dict.items():
      if count > 0:
        elements_dict_to_add[element] = count
      elif count < 0:
        elements_dict_to_subtract[element] = -count
    formula_plus=Formula(elements_dict_to_add,'None')
    formula_minus=Formula(elements_dict_to_subtract,'None')
    adduct_mass = formula_plus.get_monoisotopic_mass() - formula_minus.get_monoisotopic_mass() 
    return formula_plus, formula_minus, adduct_mass
    
    
  
  def get_multimer(self) -> int:
    return self.__multimer
  
  def get_formula_str(self) -> str:
    local_formula = f"+{str(self.__formula_plus)}" if self.__formula_plus !={} else ""
    local_formula = local_formula + f"-{str(self.__formula_minus)}" if self.__formula_minus !={} else ""
    return local_formula

  def get_formula_plus(self)-> 'Formula':
    return self.__formula_plus

  def get_formula_minus(self)-> 'Formula':
    return self.__formula_minus
  
  def get_adduct_mass(self)-> float:
    return self.__adduct_mass
  
  def get_adduct_charge(self)-> int:
    return self.__charge

  def get_adduct_charge_type(self)-> str:
    return self.__charge_type

      
    