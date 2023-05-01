import sys
sys.path.append('subformula-graphs')

import spectrum_reader
import os
from formula import Formula, FormulaGenerator
from mass_spectrum import MassSpectrum 
import matplotlib.pyplot as plt
import pandas as pd 
import constants
from timeit import timeit


#75 seconds for CHNO, 464 seconds for CHNOPF

def test_all_ms_in_directory(directory, ppm_error, default_alpha, scoring_function, output_name, delta_frag=1): #automatically tests casmi and recetox data
    
    ''' Tests, for a set of mass spectra annotated with the correct formula of the molecular ion, how well the PSG method
        retrieves the correct formula, in terms of its rank by a scoring function. Writes the output to a CSV with 4 columns,
        denoting the correct formula, rank of the correct formula, total number of candidates and the number of masses (peaks) 
        in the mass spectrum '''

    all_spectra = spectrum_reader.read_ms_from_directory(directory)
    res = []
    for index, ms_file in enumerate(all_spectra):
        metadata = ms_file.read_metadata()
        print("{} of {} files read!".format(index, len(all_spectra)))
        formula = Formula.from_string(metadata["Molecular Formula"])
        elements = set(formula.get_constituent_elements())

        test_alpha = list(default_alpha.union(elements)) #take the union of the default and the elements in the correct formula!
        test_alpha.sort(key=lambda x: constants.EXACT_MASSES.get(x))
        print(formula, test_alpha)
        spectrum = MassSpectrum(ms_file.read_mass_spectrum(), ppm_error)
        closest_mass = min(spectrum.masses, key=lambda m: abs(formula.exact_mass - m))
        
        if abs(closest_mass - formula.exact_mass) < 0.1 and len(spectrum.spectrum_dict) > 1 : #if the molecular ion is in the mass spectrum and ms has more than 1 peak:
            spectrum = spectrum.get_p_filtered_mass_spectrum(closest_mass, 0)
            
            all_annotations =  spectrum.get_spectral_annotations(test_alpha, scoring_function(spectrum), delta_frag=delta_frag)
            possible_mol_ions = [str(max(annotation,key=lambda f: f.exact_mass)) for annotation in all_annotations]
            
            try:
                rank = possible_mol_ions.index(str(formula)) + 1
            except ValueError:
                rank = None 
            
            row = (str(formula), rank, len(possible_mol_ions), len(spectrum.spectrum_dict))
            res.append(row)
            print((row))
    

    df = pd.DataFrame({"Formula": [entry[0] for entry in res], "Rank": [entry[1] for entry in res], "Candidates": [entry[2] for entry in res], "Number of Peaks": [entry[3] for entry in res]})
    df.to_csv("{}.csv".format(output_name), index=False)
    return df 



def test_all_ms_given_formula(directory, ppm, neutral_loss=False):

    ''' Tests, for a set of mass spectra annotated with the correct formula of fragment ions, how well the PSG method
        retrieves the correct fragment ion. Returns an average of the percentage of fragment ions correctly identified
        for a mass spectrum over all members of the data set'''

    all_spectra = spectrum_reader.read_ms_from_directory(directory)
    res = []
    for index, ms_file in enumerate(all_spectra):
        metadata = ms_file.read_metadata()
        ion_mode = metadata.get("Ion Mode")
        hydrogen_atom = Formula.from_string("H")
        parent_formula = Formula.from_string(metadata.get("Molecular Formula"))
        elements = set(parent_formula.get_constituent_elements())
        
        correct_annotation = ms_file.get_annotated_peaks()
        correct_formulae = [Formula.from_string(correct_annotation.get(mass)) for mass in correct_annotation]
        correct_formulae = {str(f) for f in correct_formulae if (f - hydrogen_atom).is_subformula(parent_formula)}

        spectrum = MassSpectrum(ms_file.read_mass_spectrum(), ppm)
        closest_exp_mass = min(spectrum.masses, key=lambda m: parent_formula.mass_deviation(m))

        
        if parent_formula.mass_deviation(closest_exp_mass) < 5* ppm and len(spectrum.spectrum_dict) > 1 : #if the molecular ion is in the mass spectrum and ms has more than 1 peak:
            spectrum = spectrum.get_p_filtered_mass_spectrum(closest_exp_mass, 0)
            
            if neutral_loss: #test by computing the NSG instead of the PSG
                trial_annotation = spectrum.generate_fragment_formulae_NL(parent_formula, delta_frag=1)
            else:
                trial_annotation = spectrum.generate_fragment_formulae(parent_formula, delta_frag=1)


            #Correct hydrogen counts in the generated formulae depending on ion mode
            if ion_mode == "ION_MODEPOSITIVE": 
                trial_annotation = [str(formula + hydrogen_atom) for formula in trial_annotation]
            elif ion_mode == "ION_MODENEGATIVE": 
                trial_annotation = [str(formula - hydrogen_atom) for formula in trial_annotation]
            else:
                raise Exception("Invalid ion mode!")


            num_correct = sum([formula in correct_formulae for formula in trial_annotation])

            if len(trial_annotation) != 0:
                ppv = num_correct / len(trial_annotation) #num correct / num annotations
                tpr = num_correct / len(correct_annotation) 

                print(correct_formulae)
                print([t for t in trial_annotation if t not in correct_formulae])
                print(ppv)
                res.append((ppv, tpr))
            
    
    avg_ppv = sum([e[0] for e in res])/len(res)
    avg_tpr = sum([e[1] for e in res])/len(res)
    return avg_ppv, avg_tpr



#This tests the CASMI-2016 data set and writes the output to a CSV

# test_all_ms_in_directory("mass-spectra/CASMI_2016", 1, set(["C", "H", "N", "O", "P", "F"]), lambda ms: ms.product_scoring_function, "product_chnopf", delta_frag=1)


#This tests the Recetox data set and writes the output to a CSV

# test_all_ms_in_directory("mass-spectra/Recetox", 3, set(["C", "H", "N", "O"]), lambda ms: ms.product_scoring_function, "recetox_product_chno", delta_frag=1)


#This tests the CASMI-2016 data set for the average accuracy of formula annotations for the PSG and NSG method.
#Outputs a tuple representing the average ppv and tpr

# test_all_ms_given_formula("mass-spectra/CASMI_2016", 3, neutral_loss=False)