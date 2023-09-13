#!/usr/bin/env python

"""NotFoundElement.py: Exception to raise when a chemical element is not present"""

__author__      = "Alberto Gil de la Fuente"
__copyright__   = "GPL License version 3"

class NotFoundElement(Exception):
    """
    raised when the chemical Element is not found in the periodic table
    """
    def __init__(self, input_value):
        """
        Args:
        input_value (whatever): 

        Returns:
            NotFoundElement: a new instance of a NotFoundElement Exception class
        """
        self.__input_value = input_value

    def __str__(self):
        return "The element " + self.__input_value + " does not correspond to any element"
    
    def __repr__(self):
        return "The element " + self.__input_value + " does not correspond to any element"
    
    