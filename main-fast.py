import os
import sys
import chromatogram
import formula as mf
import config
import subformula_graph
import plots


def open_file(filename):
    return chromatogram.Chromatogram(f"./{filename}")

try:
    parameters = config.parse_config("./config.json")
except KeyError:
    raise Exception("Error, either one or more parameters missing from configuration file or invalid configuration file!")

def print_formulae(formulae_list):
    for f in formulae_list:
        print("Exp Mass: {} Formula: {} Score: {} Mass Deviation: {}".format(*f))

def print_table(mat): #prints formulae list LaTeX table given the information
    print(r"\begin{table}")
    print(r"\begin{center}")
    print(r"\begin{tabular}{",end="")
    print("c"*len(mat[0]),end="")
    print(" }")
    print(r"\toprule")
    print(r"Experimental Mass & Candidate Formula & Score & Mass Deviation(ppm) \\")
    print(r"\midrule")
    for i in range(len(mat)):
        for j in range(len(mat[i])):
            if j < len(mat[i]) - 1:

                if j == 1:
                    print(r"\ce{" + "{}".format(mat[i][j]),end="" + r"} & ")
                else:
                    print("{} &".format(mat[i][j]),end=" ")
            else:
                print("{} \\\ ".format(mat[i][j]))
    print(r"\bottomrule")
    print(r"\bottomrule")
    print(r"\end{tabular}")
    print(r"\end{center}")
    print(r"\end{table}")


data = open_file("WF12011.CDF")

plots.plot_chromatogram(data, parameters.peak_sensitivity)
n = 5 # int(input("Enter Scan Number: "))
ms = data.get_raw_MS(n).binned_MS().conventional_norm().intensity_cutoff(parameters.intensity_factor, 999)
plots.plot_mass_spectrum(ms)   
ions = subformula_graph.find_molecular_ion(ms, *parameters.mass_range, parameters.alphabet, parameters.ppm_error, parameters.heuristic)
print("\n\nResults:\n")
print_formulae(ions)
