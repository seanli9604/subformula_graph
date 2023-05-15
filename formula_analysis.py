import sys
sys.path.append('subformula-graphs')

import spectrum_reader
import mass_spectrum
from mass_spectrum import MassSpectrum
import spectrum_reader
from visualisation import MSVisualiser
import matplotlib.pyplot as plt  
import pandas as pd
from timeit import timeit

#The alphabet for the semiochemical data set
alpha = ["C","H","N","O","S"]

df = pd.read_csv("mass-spectra/semiochemical-dataset.csv") #reads semiochemical data set from csv

#The semiochemical dataset, denoted by a list of tuples representing the compound ID filename, the scan number and the 
#correct molecular formula as identified by an expert analyst
data = [(row[0], row[1], row[2], row[3]) for row in zip(df['Compound ID'], df['File Name'], df['Scan Number'], df['Correct Formula'])] #list of rows of the data set 

#time for analysis: 44.55813224599842

#Change this to change the scoring function used!
func = lambda ms: ms.vertex_scoring_function


def get_rank_one_case(row, scoring_function = func, ppm_error=10, intensity_bound=0.1, lower_bound=130, upper_bound=250): #returns the rank of the correct annotation
    compound_id, filename, scan_number, correct_annotation = row
    ms_dict = spectrum_reader.read_ms_from_file("mass-spectra/chromatograms/{}.CDF".format(filename)).get_baseline_subtracted_mass_spectrum(scan_number)
    ms = MassSpectrum(ms_dict, ppm_error).conventional_norm()
    scoring_func = func(ms)
    annotations = ms.compute_most_likely_molecular_ion(alpha, scoring_func, lower_bound, upper_bound, cutoff_val=intensity_bound)
    candidate_molecular_ions = [str(a[-1]) for a in annotations]
    

    if correct_annotation in candidate_molecular_ions:
        rank = candidate_molecular_ions.index(correct_annotation) + 1
        print(ms.get_formula_dict(annotations[rank-1]))
    else:
        rank = 0

    return rank
    

def get_rank_all_cases(dataset): #iterate through the whole dataset and returns a list of ranks 
    res = []
    for row in dataset:
        rank = get_rank_one_case(row)
        print(rank)
        res.append(rank)
    return res 


def make_2D_fragment_plot(row, scoring_function = func, ppm_error=10, intensity_bound=0.1, lower_bound=130, upper_bound=250, num_plots=5): #generates a 2DFP for the given entry
    compound_id, filename, scan_number, correct_annotation = row
    ms_dict = spectrum_reader.read_ms_from_file("mass-spectra/chromatograms/{}.CDF".format(filename)).get_baseline_subtracted_mass_spectrum(scan_number)
    ms = MassSpectrum(ms_dict, ppm_error).conventional_norm()
    scoring_func = func(ms)
    annotations = ms.compute_most_likely_molecular_ion(alpha, scoring_func, lower_bound, upper_bound, cutoff_val=intensity_bound)

    for index, annotation in enumerate(annotations[:num_plots]):
        candidate_formula = max(annotation, key=lambda f: f.exact_mass)
        vis = MSVisualiser(ms, annotation)
        vis.plot_ms_as_2D_fragment()
        #plt.savefig(new_dir + "/" + "{}_rank_{}_correct_formula_{}.png".format(candidate_formula, index+1, correct_formula),dpi=96)
        plt.show()

    return


def df_to_latex_table(df, formulae=[], captiontext="Insert Caption!"): #converts a pandas dataframe to a printed latex table 

    names = [name for name in df]

    print(r"\begin{table}")
    print(r"\begin{center}")
    print(r"\begin{tabular}{",end="")
    print("c"*len(names),end="")
    print(" }")
    print(r"\hline") #toprule
    print(" & ".join(names) + r"  \\" )
    print(r"\hline") #midrule

    for index, row in df.iterrows():
        formatted_row = [str(row[name]) if index not in formulae else "\ce{" + row[name] + "}" for index, name in enumerate(names)]
        newline = " & ".join(formatted_row) + r"  \\"
        print(newline)

    print(r"\hline") #bottomrule
    print(r"\hline") #bottom rule
    print(r"\end{tabular}")
    print(r"\caption{" + captiontext + r"}") 
    print(r"\end{center}")
    print(r"\end{table}")
    print("\n\n")


def analyse_one_case(row, scoring_function = func, ppm_error=10, intensity_bound=0.1, lower_bound=130, upper_bound=250):

    ''' Generates a pandas dataframe, and some caption text containing information about each formula annotation'''

    compound_id, filename, scan_number, correct_annotation = row
    ms_dict = spectrum_reader.read_ms_from_file("mass-spectra/chromatograms/{}.CDF".format(filename)).get_baseline_subtracted_mass_spectrum(scan_number)
    ms = MassSpectrum(ms_dict, ppm_error).conventional_norm()
    scoring_func = func(ms)
    annotations = ms.compute_most_likely_molecular_ion(alpha, scoring_func, lower_bound, upper_bound, cutoff_val=intensity_bound)
    tic = sum(ms.intensities) #total ion current 

    candidate_molecular_ions = [str(a[-1]) for a in annotations] #last item of annotation list = parent ion
    formula_masses = [a[-1].exact_mass for a in annotations]
    scores = [scoring_func(a) for a in annotations]

    all_formula_dicts = [ms.get_formula_dict(a) for a in annotations]
    percent_spectrum_explained = [sum([formula_dict.get(f) for f in formula_dict]) * 100 / tic for formula_dict in all_formula_dicts]
    
    if correct_annotation in candidate_molecular_ions:
        rank = candidate_molecular_ions.index(correct_annotation) + 1

    else:
        rank = 0

    num_rows = max([5, rank])
    df = pd.DataFrame({"Rank": [i+1 for i in range(num_rows)], "Formula": [i for i in candidate_molecular_ions[:num_rows]],
    "Monoisotopic Mass": [round(i,3) for i in formula_masses[:num_rows]], "Score": [round(i,3) for i in scores[:num_rows]], 
    "Spectrum Explained (\\%)": [round(i, 2) for i in percent_spectrum_explained[:num_rows]]})
    captiontext = "A list of candidate molecular formulae annotated to the mass spectrum of Compound {}. The correct molecular formula is \\ce{{{}}}, and its rank is {} using the product score metric.".format(compound_id, correct_annotation, rank)    
    return df, captiontext 


def generate_all_latex_tables(dataset):
    
    ''' Given the dataset, generate LaTeX tables (with captions) summarising 
         the annotations generated by the program '''
    
    for row in dataset:
        df, caption = analyse_one_case(row)
        df_to_latex_table(df, formulae=[1], captiontext=caption)
    

get_rank_all_cases(data)