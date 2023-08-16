import re

def atom_dict_from_str(formula_str: str) -> 'Dict':

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
      
    return elements
