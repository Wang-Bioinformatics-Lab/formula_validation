#!/usr/bin/env python

"""
This Python module contains not only the class FormulaPrediction, but also the test of this FormulaPrediction class.

@contents :  This Python module contains not only the class FormulaPrediction, but also the test of this FormulaPrediction class. 
It uses the R function from https://github.com/mjhelf/MassTools/blob/054ae740fa63aa16d6757bc6080b0ec380680a76/R/Functions_MF_prediction.R to check if a monoisotopic mass can be explainable from a molecular formula
@project :  N/A
@program :  N/A
@file :  FormulaPrediction.py
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
from Formula import Formula

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

robjects.r('''
    install.packages("devtools")
    devtools::install_github("mjhelf/MassTools")

    library(MassTools)
''')

mz = 200
charge = 1
ppm = 20
min_carbons = 8
max_carbons = 9
min_hydrogens = 4
max_hydrogens = 9
top = 3

min_elements = f"C{min_carbons}H{min_hydrogens}"
max_elements = f"C{max_carbons}H{max_hydrogens}"

elements_list=["C", "H", "N", "O", "S"]
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
robjects.r('mz <- {}'.format(mz))
robjects.r('z <- {}'.format(charge))
robjects.r('ppm <- {}'.format(ppm))
robjects.r('top <- {}'.format(top))
robjects.r('mfs <- calcMF(mz = mz, z = z, ppm = ppm, top = top, elements = elements, Filters = {})'.format(filter_r_code))
mfs = robjects.r('mfs')
