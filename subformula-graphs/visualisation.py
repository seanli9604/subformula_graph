'''
Author: Sean Li
Date: May 1, 2023

This file generates the 2-dimensional fragment plots from a particular set of whole-spectrum 
annotation.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any 
later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. 
If not, see <https://www.gnu.org/licenses/>.
'''

from mass_spectrum import MassSpectrum
from formula import Formula
import seaborn as sns
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd


def cluster(l, clustersize): #cluster continguous elements of a list e.g [1,2,4,5] -> [[1,2], [4,5]] if clustersize = 1
    res = [[l[0]]]
    for i in l[1:]:
        if i < res[-1][-1] + clustersize:
            res[-1].append(i)
        else:
            res.append([i])
    return res


def cluster_formulae(formula_list): #cluster formulae which differ by an integer number of masses 
    formula_list = sorted(formula_list, key=lambda f: f.exact_mass) #sort list in order of mass first
    res = [[formula_list[0]]]
    for f in formula_list[1:]:
        if f.exact_mass < res[-1][-1].exact_mass + 1.1:
            res[-1].append(f)
        else:
            res.append([f])
    return res



class MSVisualiser:


    def __init__(self, input_spectrum, annotation):
        self.num_peaks = len(input_spectrum.masses) 
        input_spectrum = input_spectrum.conventional_norm()
        
        self.mass_spectrum = input_spectrum
        self.spectrum_dict = input_spectrum.spectrum_dict
        
        annotation.sort(key=lambda f: f.exact_mass, reverse=True)
        self.annotation = annotation
        self.tot_intensity = sum([self.spectrum_dict.get(m) for m in self.spectrum_dict])
        self.formula_dict = input_spectrum.get_formula_dict(annotation)
    
    
    def get_serialised_ms(self, compression=0.1): #"serialise" a mass spectrum. That is, represent a higher abundance as a larger count in the number of masses
        arraylist = []
        for mass in self.spectrum_dict:
            numcounts = int(self.spectrum_dict[mass] / compression)
            counts = np.full((numcounts,), int(round(mass)))
            arraylist.append(counts)
        return np.concatenate(arraylist)
    

    def generate_formula_labels(self): #generate all formula labels from the annotations by removing ions differing by just a hydrogen
        clustered_formulae = cluster_formulae(self.annotation)
        return [max(cluster, key=lambda f: self.formula_dict.get(str(f))) for cluster in clustered_formulae]


    def stagger_labels(self, clustersize): #Returns a dictionary mapping masses to a padding value such that the labels are staggered
        formula_masses_dict = {formula.nom_mass: formula for formula in self.generate_formula_labels()}
        clusters = cluster([f for f in formula_masses_dict], clustersize)
        stagger_dict = {}
        for group in clusters:
            for index, mass in enumerate(group):
                stagger_dict[mass] = 1 + 10 * (index)
        return stagger_dict


    def plot_ms_as_2D_fragment(self, clustersize=5): 

        ''' Plot the mass spectrum as a 2D fragment plot (2DFP) '''
        
        formula_masses_dict = {formula.nom_mass: formula for formula in self.annotation}
        PSG_edges = self.mass_spectrum.make_edges_for_annotation(self.annotation) #edges of the PSG
        stagger_dict = self.stagger_labels(clustersize)

        ms = self.get_serialised_ms()
        ms_max = int(self.annotation[0].exact_mass) + 6
        mass_range = ms_max - min(ms)

        df = pd.DataFrame({"index": np.arange(len(ms)), "m/z": ms})
        ms_filtered = np.array([m for m in ms if m in formula_masses_dict])
        df_filtered = pd.DataFrame({"index": np.arange(len(ms_filtered)), " ": ms_filtered})

        spectrum_bins = np.linspace(min(ms), ms_max + 6, mass_range*2 + 12, endpoint=True) #add some padding to the end of the histogram so that the formulae will fit!
        grid = sns.JointGrid(height=10, ratio=4, space=.25, marginal_ticks=False)
        subplot_x = sns.histplot(x=df['m/z'], fill=True, linewidth=0, ax=grid.ax_marg_x, bins= spectrum_bins, color = "red")
        subplot_y = sns.histplot(y=df['m/z'], fill=True, linewidth=0, ax=grid.ax_marg_y, bins = spectrum_bins, color = "red")
        formula_plot = sns.scatterplot(x=df_filtered[' '], y=df_filtered[' '], ec="blue", fc="blue", s=25, linewidth=0, ax=grid.ax_joint)  
        
        #add the green dots representing the PSG
        subformula_df = pd.DataFrame({"formula": [edge[1].nom_mass for edge in PSG_edges], "subformula": [edge[0].nom_mass for edge in PSG_edges]})
        sns.scatterplot(x=subformula_df['subformula'], y=subformula_df['formula'], ec="green", fc="green", s=25, linewidth=0, ax=grid.ax_joint)

        for i in range(mass_range*2): #set peaks matching annotations to blue
            mass = min(ms) + i / 2            

            if int(round(mass + 0.001, 0)) in formula_masses_dict:
                grid.ax_marg_x.patches[i].set_facecolor("blue")
                grid.ax_marg_y.patches[i].set_facecolor("blue")
        

        refined_formula_masses_dict = {formula.nom_mass: formula for formula in self.generate_formula_labels()}
        for mass in refined_formula_masses_dict: #add formula labels
            grid.refline(x=mass, y=mass, marginal=False, color='.5', alpha=0.2)
            grid.ax_joint.annotate(text = refined_formula_masses_dict[mass].latex_formula(), xy = (mass + 1 + stagger_dict[mass], mass -1 ), fontsize=10, fontweight="bold") #+ stagger_dict[mass]

        
        product_score = round(self.mass_spectrum.product_scoring_function(self.annotation), 3)
        largest_fragment = self.annotation[0] - self.annotation[1]
        frag_string =  str(largest_fragment.nom_mass) + "({})".format(str(largest_fragment))
        percentage_spec_explained = round(sum([self.formula_dict.get(str(f)) for f in self.annotation]) * 100 / self.tot_intensity, 1)
        grid.ax_joint.annotate(text=" Score: {} \n M1 - M2: {} \n % Explained: {} ".format(product_score, frag_string, percentage_spec_explained), xy = (ms_max*0.8, min(ms)*1.2), fontsize=14, fontweight="bold", color = "green")


        
    


    


    

    
