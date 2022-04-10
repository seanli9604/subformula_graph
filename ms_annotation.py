import os
import sys
import netCDF4 as nc
import formula as mf
import frag_tree as ft
import netcdf
import networkx as nx
import matplotlib.pyplot as plt 


def file_loop():
    while True:
        try:
            a = input("Please enter the name of your data file (case sensitive), making sure to include file extensions (e.g MyFile.cdf): ")
            ### need to add ways to handle other file types ###
            data = netcdf.Chromatogram(os.path.dirname(sys.argv[0]) + "\\" + a)
            return data
    
        except FileNotFoundError:
            print("File not found! Make sure you have included the file extension in the name and the data file is in this folder! \n")



def plot_chromatogram(data:netcdf.Chromatogram, intensity_cutoff):
    y = data.ion_current
    t = [data.times[i]/60 for i in range(data.maxscan)]
    plt.plot(t,y)

    peaks = data.get_TIC_peaks(intensity_cutoff)

    for p in peaks:
        plt.annotate("{} \nRT={}".format(p, round(t[p],2)), # this is the text
                (t[p], y[p]), # these are the coordinates to position the label
                fontsize=8,
                textcoords="offset points", # how to position the text
                xytext=(0,2), # distance from text to points (x,y)
                ha='center') # horizontal alignment can be left, right or center

    plt.ylabel("Total ion current")
    plt.xlabel("Time (minutes)")
    plt.show()


def plot_mass_spectrum(spectrum):
    max_mass = int(spectrum[-1][0]) + 50
    d = {int(round(p[0],0)): p[1] for p in spectrum}
    mz = [i for i in range(max_mass)]
    a = [d[i] if i in d else 0 for i in range(max_mass)]
    plt.bar(mz, a, color ='blue',
        width = 0.4)
    
    for p in spectrum:
        plt.annotate(p[0], # this is the text
                (int(round(p[0],0)), p[1]),
                fontsize=8, 
                textcoords="offset points", # how to position the text
                xytext=(0,2), # distance from text to points (x,y)
                ha='center') # horizontal alignment can be left, right or center

    plt.xlabel("m/z")
    plt.ylabel("Abundance")
    plt.title("Mass spectrum")
    plt.show()


def plot_annotated_mass_spectrum(spectrum, graph):
    max_mass = int(spectrum[-1][0]) + 50
    d = {int(round(p[0],0)): p[1] for p in spectrum}
    mz = [i for i in range(max_mass)]
    a = [d[i] if i in d else 0 for i in range(max_mass)]
    plt.bar(mz, a, color ='blue',
        width = 0.4)
    
    fs = ft.return_formulae(ft.extract_mass(spectrum),graph,to_dict=True)
    for p in spectrum:
        if int(round(p[0],0)) in fs:
            plt.annotate(fs[int(round(p[0],0))], # this is the text
                (int(round(p[0],0)), p[1]),
                fontsize=8, 
                textcoords="offset points", # how to position the text
                xytext=(0,2), # distance from text to points (x,y)
                ha='center') # horizontal alignment can be left, right or center

    plt.xlabel("m/z")
    plt.ylabel("Abundance")
    plt.title("Mass spectrum")
    plt.show()


def MS_loop(data:netcdf.NetCDF4ScanData, max_mass:int, noise_factor:float, int_cutoff:float):
    while True:
        try:
            print("The maximum scan number is {}".format(data.maxscan))
            n = input("Please enter a scan number from which you want to retrieve your mass spectrum: ")
            n = int(n.strip())
            try: 
                print("Scan number received! Plotting mass spectrum..")
                background = data.background(max_mass)
                ms = data.make_raw_MS(n).de_noise(background,noise_factor,int_cutoff)
                plot_mass_spectrum(ms.get_MS())
                answer = input("Do you want to select this mass spectrum for analysis? (y/n): ")
                if answer == "y" or answer == "Y" or answer == "yes" or answer == "YES":
                    return n

            except ValueError:
                print("Error, there doesn't seem to be a peak here!")
        except ValueError:
            print("Error, need to input integer!")
        
        

def formula_loop(graphs, spec):
    while True:
        try:
            print("There are {} sets of annotations possible (1 - {})".format(len(graphs),len(graphs)))
            num = int(input("Please select which set you want to view: "))
            l = ft.return_formulae(spec, graphs[num -1])
            
            for formula in l:
                print("Mass: {}                  Formula: {}".format(formula[1], formula[0]))
            
            answer = input("Do you want to choose this set of annotations? (y/n): ")
            if answer == "y" or answer == "Y" or answer == "yes" or answer == "YES":
                return num -1



        except ValueError or IndexError:
            print("Invalid input! Please enter a number between 1 and {}".format(len(graphs)))

############################# main loop #############################

print("""Welcome to the mass spectrometry visualisation and formula annotation program!
This program will walk you through your MS spectral analysis. \n""")\

data = file_loop()

print("Data successfully read!")
while True:
    print("Visualising chromatogram...")

    plot_chromatogram(data,10000)

    n = MS_loop(data,500,3,0.005)
    print("You have chosen to analyse scan no {}".format(n))
    print("Retrieving spectrum...")
    
    background = data.background(500)
    ms = data.make_raw_MS(n).de_noise(background,3,0.005)

    print("Computing formulae...")
    g = ft.analyse_MS(ms.masses)

    print("Formula computation complete!")
    num  = formula_loop(g, ms.masses)

    ans = input("Do you want to see the subformula graph of this annotation? (y/n): ")
    if ans == "y" or ans == "Y" or ans == "yes" or ans == "YES":
        ft.draw_graph(g[num][0])

    print("Displaying annotated mass spectrum...")
    plot_annotated_mass_spectrum(ms.get_MS(),g[num])
    print("Analysis complete!")
    ans = input("Do you want to do another analysis? (y/n): ")
    if ans == "y" or ans == "Y" or ans == "yes" or ans == "YES":
        pass
    else:
        break




