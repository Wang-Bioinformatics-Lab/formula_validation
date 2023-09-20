#!/usr/bin/env python

"""
This Python module contains the tests of the class Adduct. 

@contents :  This Python module contains the tests of the class Adduct. 
@project :  N/A
@program :  N/A
@file :  test_adduct.py
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
from formula_validation.Adduct import Adduct
from formula_validation.Element import Element_type, element_weights
from formula_validation.IncorrectFormula import IncorrectFormula
from formula_validation.IncorrectAdduct import IncorrectAdduct
from formula_validation.NotFoundElement import NotFoundElement



class TestAdduct(unittest.TestCase):
  def setUp(self):
    # Set up any common resources needed for tests
    pass

  def tearDown(self):
    # Clean up after each test
    pass

  def test_adduct_m_plus_h(self):
    adduct_str1 = '[M+H]+'
    adduct_1 =Adduct(adduct_str1)
    self.assertEqual(adduct_1.get_adduct_charge(), 1)
    self.assertEqual(adduct_1.get_adduct_charge_type(), '+')
    self.assertEqual(adduct_1.get_multimer(), 1)
    self.assertAlmostEqual(adduct_1.get_adduct_mass(), 1.0078, delta=1e-3)


  def test_adduct_2m_plus_whatever(self):
    adduct_str1 = '[2M+HCOOH-H]-'
    adduct_1 =Adduct(adduct_str1)
    self.assertEqual(adduct_1.get_adduct_charge(), 1)
    self.assertEqual(adduct_1.get_adduct_charge_type(), '-')
    self.assertEqual(adduct_1.get_multimer(), 2)
    self.assertAlmostEqual(adduct_1.get_adduct_mass(), 44.9976, delta=1e-3)

  def test_adduct_5m_plus_whatever(self):
    adduct_str1 = '[5M+Ca]2+'
    adduct_1 =Adduct(adduct_str1)
    self.assertEqual(adduct_1.get_adduct_charge(), 2)
    self.assertEqual(adduct_1.get_adduct_charge_type(), '+')
    self.assertEqual(adduct_1.get_multimer(), 5)
    self.assertAlmostEqual(adduct_1.get_adduct_mass(), 39.9625, delta=1e-3)

  def test_adduct_m_plus_3H2O(self):
    adduct_str1 = '[M-3H2O+2H]2+'
    adduct_1 =Adduct(adduct_str1)
    self.assertEqual(adduct_1.get_adduct_charge(), 2)
    self.assertEqual(adduct_1.get_adduct_charge_type(), '+')
    self.assertEqual(adduct_1.get_multimer(), 1)
    self.assertAlmostEqual(adduct_1.get_adduct_mass(), -52.016, delta=1e-3)

  def test_adduct_m_minus2HPlusK(self):
    adduct_str1 = '[M-2H+K]-'
    adduct_1 =Adduct(adduct_str1)
    self.assertEqual(adduct_1.get_adduct_charge(), 1)
    self.assertEqual(adduct_1.get_adduct_charge_type(), '-')
    self.assertEqual(adduct_1.get_multimer(), 1)
    self.assertAlmostEqual(adduct_1.get_adduct_mass(), 36.94806, delta=1e-3)

  def test_adduct_m_minus2HPlusK_neutral(self):
    adduct_str1 = '[M-2H+K]'
    adduct_1 =Adduct(adduct_str1)
    self.assertEqual(adduct_1.get_adduct_charge(), 0)
    self.assertEqual(adduct_1.get_adduct_charge_type(), '')
    self.assertEqual(adduct_1.get_multimer(), 1)
    self.assertAlmostEqual(adduct_1.get_adduct_mass(), 36.94806, delta=1e-3)


   

if __name__ == "__main__":
    unittest.main()