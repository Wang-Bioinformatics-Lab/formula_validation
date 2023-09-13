#!/usr/bin/env python

"""NotFoundElement.py: Exception to raise when a chemical element is not present"""

__author__      = "Alberto Gil de la Fuente"
__copyright__   = "GPL License version 3"

class IncorrectAdduct(Exception):
    """
    raised when the adduct format is not valid. An adduct is represented by '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+'
    """
    def __init__(self, input_value):
        """
        Args:
        input_value (whatever): 

        Returns:
            IncorrectAdduct: a new instance of a IncorrectAdduct Exception class
        """
        self.__input_value = input_value

    def __str__(self):
        return "The adduct " + self.__input_value + " does not correspond to a correct adduct"
    
    def __repr__(self):
        return "The adduct " + self.__input_value + " does not correspond to a correct adduct"
    
    