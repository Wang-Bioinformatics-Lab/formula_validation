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
# libraries to make requests to chemcalc
import urllib3
import json

from Element import Element_type, element_weights
from IncorrectFormula import IncorrectFormula
from NotFoundElement import NotFoundElement
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.rinterface_lib.sexp import NULLType

class Formula:
  __electron_weight=0.00054858

  __connection_to_chemcalc_ws = urllib3.PoolManager()
  __ccurl = 'https://www.chemcalc.org/chemcalc/mf' 
  
  __default_ppm = 50

  robjects.r('''
    install.packages("devtools")
    devtools::install_github("mjhelf/MassTools")

    library(MassTools)
  ''')

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
        self.__adduct = adduct
    else: 
      raise IncorrectFormula("The adduct " + str(adduct) + " is not a valid adduct with the format: '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+' ")
    
    self.__monoisotopic_mass_with_adduct = self.__calculate_monoisotopic_mass_with_adduct()

    
  def __eq__(self, other):
    if not isinstance(other, Formula):
      return False
    return self.__elements == other.get_elements()

  def __str__(self):
    result_string = "".join(str(key.name) + str(value) for key,value in self.__elements.items())
    return result_string
  
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
      return Formula(new_formula_dict, self.__adduct)
      
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
      return Formula(new_formula_dict, self.__adduct)
    else:
        raise IncorrectFormula("other should be a formula and is a " + str(type(Formula)))
    

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
    if re.search(r'(?<![A-Z])[a-z]', formula_str):
      raise IncorrectFormula("The formula contains elements that are not chemical Elements")
    # check for any character that is not a capital letter, a lowercase letter, or a number
    if re.search(r'[^a-zA-Z0-9]', formula_str):
      
      raise IncorrectFormula("The formula contains parenthesis or brackets")

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
      simple_formula = Formula.formula_from_str_hill(formula_str, adduct)
      return simple_formula
    except Exception as e:
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
    

  def check_monoisotopic_mass(self, external_mass: Union[float,int], mass_tolerance_in_ppm: Union[int, float] = __default_ppm) -> bool:
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
  
  def check_monoisotopic_mass_with_adduct(self, external_mass: Union[float,int], mass_tolerance_in_ppm: Union[int, float] = __default_ppm) -> bool:
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


  def ppm_to_absolute(reference_monoisotopic_mass: Union[float,int], ppm: Union[float,int] = __default_ppm) -> float:
    """
    Args:
      reference_monoisotopic_mass (numeric): monoisotopic mass of reference
      ppm (numeric): ppm of the reference monoisotopic mass
    Returns: 
      the absolute value of the ppm calculated
    Raise: 
      a exception if reference_monoisotopic_mass or ppm are not numbers
    """
    return (reference_monoisotopic_mass / 1000000.0) * ppm
  
  def check_possible_fragment_mz(self, fragment_mz: Union[float, int], ppm: Union[float, int] = __default_ppm):
    """
    Args:
      self (Formula): formula to check if the fragment could be generated by this formula according to the substructure and the adduct. If the adduct has a number of charges > 1, it will check the adducts for double and single charge
      fragment_mz (numeric): mz_mass to check if it can come from the structure
      ppm (numeric): ppm of the reference monoisotopic mass to check the potential formulas
    Returns: 
      a boolean specifying if the fragment m_z can be explained from the formula and the adduct
    Raise: 
      a exception if fragment_mz or ppm are not numbers
    """

    mz = 200
    charge = 1
    
    top = 5
    final_formula = self
    if self.__adduct is not None:
      charge= self.__adduct.get_adduct_charge()
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
    robjects.r('z <- {}'.format(charge))
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
    Returns: 
      the percentage of intensity of fragments explained according to the formula and adduct
    Raise: 
      a exception if fragment_mz or ppm are not numbers
    """
    
    total_fragment_intensities = sum(value for value in fragments_mz_intensities.values() if isinstance(value, (int, float)))
    explained_intensities = 0
    for mz, intensity in fragments_mz_intensities.items():
      if self.check_possible_fragment_mz(mz,ppm):
        explained_intensities +=intensity
    return explained_intensities / total_fragment_intensities