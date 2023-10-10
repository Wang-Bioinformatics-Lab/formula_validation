from formula_validation.Formula import Formula
import csv


def main():
    adducts_file = 'adducts_maricruz.txt'
    updated_data = []
    adducts = []
    column_monoisotopic_mass_fixed_name = 'monoisotopic_mass_fixed'

    # Open the text file using a context manager
    with open(adducts_file, mode='r') as adducts_file:
        # Read each line from the file
        for adduct in adducts_file:
            adducts.append(adduct.strip())
    
    formulas_file = 'tabla_formulas.csv'
    formulas_file_corrected = 'tabla_formulas_corrected.csv'

    with open(formulas_file, mode='r', newline='') as formulas_file, \
        open(formulas_file_corrected, mode='w', newline='') as new_file:
        csv_reader = csv.DictReader(formulas_file, delimiter='|')
        fieldnames = csv_reader.fieldnames
        fieldnames.append(column_monoisotopic_mass_fixed_name)
        fieldnames = fieldnames + adducts
        csv_writer = csv.DictWriter(new_file, fieldnames=fieldnames, delimiter='|')
        csv_writer.writeheader()

        
        for row in csv_reader:
            formula_str = row['FORMULA']
            formula = Formula.formula_from_str(formula_str,None)
            monoisotopic_weight = formula.get_monoisotopic_mass()
            row[column_monoisotopic_mass_fixed_name] = monoisotopic_weight
            for adduct in adducts:
                formula_adduct = Formula.formula_from_str(formula_str,adduct)
                monoisotopic_weight_adduct = formula_adduct.get_monoisotopic_mass_with_adduct()
                row[adduct] = monoisotopic_weight_adduct
            
            csv_writer.writerow(row)
        
            
if __name__ == "__main__":
    main()