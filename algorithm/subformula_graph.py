import itertools
import matplotlib.pyplot as plt 
import networkx as nx
import algorithm.chromatogram as chromatogram 

import algorithm.formula as mf 
from typing import List 

Vector = List[float]

alpha = ["C","H","N","O", "S"] 

def draw_graph(graph:nx.digraph):
    nx.draw_networkx_nodes(graph, pos=nx.circular_layout(graph), node_size=1000, node_color="#00ff00")
    nx.draw_networkx_labels(graph, pos=nx.circular_layout(graph), font_size=6)
    nx.draw_networkx_edges(graph, pos=nx.circular_layout(graph))
    #nx.draw_circular(graph,with_labels=True, node_size=1000)
    plt.show()
    return

class ExpMass:
    
    def __init__(self, mass, error): #error in ppm!
        self.mass = mass
        self.rounded_mass = int(round(mass,0))
        self.ppm_error = error
        self.abs_error = error*mass/10**6
        
    def __repr__(self):
        return "{} ± {}".format(self.mass, self.abs_error)
    
    def possible_formulae(self, alphabet:List[str], heuristics=None, parent=None)->List[mf.Formula]:
        e = mf.MFGenerator(alphabet, self.rounded_mass)
        
        if parent != None:
            e.parental_restriction(parent)
        elif heuristics != None:
            e.custom_heuristic(heuristics)
        else:
            e.default_heuristic()

        return  e.get_formula_list(self.mass, self.ppm_error, DBE_restriction = lambda x: x.dbe() >= 0)


class CandidateGenerator:

    def __init__(self, molecular_ion_mass:float, masses:List[float], alphabet:List[str], ppm_error:float):
        self.alphabet = alphabet
        self.ppm_error = ppm_error
        self.masses = [ExpMass(mass, self.ppm_error) for mass in masses]
        self.molecular_ion_mass = ExpMass(molecular_ion_mass, 3*self.ppm_error) # putative mol formula -- 3 standard deviations

    def read_in_spectrum(molecular_ion_mass: int, MS:chromatogram.MassSpectrum, alphabet:List[str], ppm_error:int):
        return CandidateGenerator(molecular_ion_mass, MS.masses, alphabet, ppm_error)
    
    def manual_selection(self, heuristics=None):
        candidates = [formula for formula in self.molecular_ion_mass.possible_formulae(self.alphabet, heuristics) if formula.dbe().is_integer()]
        res = []

        for c in candidates:
            frag_set, assigned_masses = [],[]
            for mass in self.masses:
                if mass.mass < self.molecular_ion_mass.mass:
                    formula_assignments = mass.possible_formulae(self.alphabet, parent=c)
                    if len(formula_assignments) != 0:
                        frag_set.append(formula_assignments)
                        assigned_masses.append(mass)
            res.append(CandidateSet(c, frag_set, assigned_masses + [self.molecular_ion_mass]))

        return sorted(res, key=lambda x: x.score_candidate(), reverse=True)


class CandidateSet: #set of formula candidates
    
    def __init__(self, molecular_ion:mf.Formula, fragments:List[List[mf.Formula]], masses:ExpMass):
        self.molecular_ion = molecular_ion
        self.fragments = fragments
        self.masses = masses
        self.first_fragment = [i[0] for i in self.fragments]
    
    def __repr__(self):
        pair = zip(self.first_fragment + [self.molecular_ion], self.masses)
        exp_error = [p[0].mass_deviation(p[1].mass) for p in pair]
        table = ""
        for exp_mass, formula, error in zip(self.masses, self.first_fragment + [self.molecular_ion] , exp_error):
            m = round(float(exp_mass.mass), 4)
            err = round(float(error), 3)
            table += "Experimental Mass: {}       Formula: {}        Deviation(ppm): {}\n".format(m, formula, err)
        return table
    
    def get_subformula_graph(self, complement=False):
        G = nx.DiGraph()
        ion_list = self.first_fragment + [self.molecular_ion]
        G.add_nodes_from(ion_list)
        if complement:
            G.add_edges_from([(ion_list[i], ion_list[j]) for i in range(len(ion_list)) for j in range(i) if not mf.Formula.is_subformula(ion_list[i],ion_list[j])])
        else:
            G.add_edges_from([(i,j) for i in ion_list for j in ion_list if mf.Formula.is_subformula(i,j) and not mf.Formula.same_formula(i,j)])       
        return G
    
    def score_candidate(self):
        G = self.get_subformula_graph(complement=False)
        V = nx.number_of_nodes(G)
        return nx.number_of_edges(G)
    
    def mass_score(self):
        return sum([t[0].mass_deviation(t[1].mass) for t in zip(self.first_fragment, self.masses)])/len(self.masses)
    
    def optimise_fragment_formulae(self):

        def wrap(l):
            return [[i] for i in l]
         
        all_fragment_combinations = [ wrap(frag) for frag in itertools.product(*self.fragments)]
        candidate_list = [CandidateSet(self.molecular_ion, frag, self.masses) for frag in all_fragment_combinations]
        return max(candidate_list, key=lambda x: x.score_candidate())
        

def find_molecular_ion(MS, lower, upper, alpha, ppm_error, heuristics=None):
    res = []
    candidate_masses = MS.find_candidate_masses(lower, upper)
    print("Found {} possible molecular masses in given mass range!".format(len(candidate_masses)))
    for ind, m in enumerate(candidate_masses):
        print("Scoring {} out of {} masses!".format(ind+1, len(candidate_masses)))
        candidates = CandidateGenerator.read_in_spectrum(m, MS, alpha, ppm_error).manual_selection(heuristics)
        for c in candidates:
            mol = c.molecular_ion
            rounded_mass = round(float(m),4)
            res.append((rounded_mass, mol, c.score_candidate(), round(mol.mass_deviation(m),2)))
    return sorted(res, key=lambda x: x[2], reverse=True)




