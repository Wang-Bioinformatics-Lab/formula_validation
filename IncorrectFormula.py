#!/usr/bin/env python

"""IncorrectFormula.py: Exception to raise when a chemical formula is not correct"""

__author__      = "Alberto Gil de la Fuente"
__copyright__   = "GPL License version 3"

class IncorrectFormula(Exception):
    """
    raised when the chemical Formula is not valid. A chemical formula is defined by a set of elements in the sequence ([CHEMICAL_ELEMENT][NUMBER_OCCURRENCES])+
    """
    def __init__(self, input_value):
        """
        Args:
        input_value (whatever): 

        Returns:
            IncorrectFormula: a new instance of a IncorrectFormula Exception class
        """
        self.__input_value = input_value

    def __str__(self):
        return "The Formula " + self.__input_value + " does not correspond to a correct formula"
    
    def __repr__(self):
        return "The Formula " + self.__input_value + " does not correspond to a correct formula"
    
    