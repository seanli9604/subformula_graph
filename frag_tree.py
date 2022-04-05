import functools
import numpy as np
import math 
import matplotlib.pyplot as plt 
import networkx as nx
import netcdf as nc
import formula as mf 
from typing import List 

Vector = List[float]

alpha = ["C","H","N","O","P","S"]

############################################set tolerances############################################ 

def edge_tol(m1:float, m2:float, ppm = 10)->float:
    std1, std2 = m1*ppm/10**6, m2*ppm/10**6
    return 3*math.sqrt(std1**2 + std2**2) #3 sigmas 


def vertex_tol(m:float, ppm=20)->float:
    return 3*m*ppm/10**6

############################################Spectrum handling############################################ 

def extract_mass(MS:List[tuple])->Vector: #from a mass spectrum return a list of masses
    return [i[0] for i in MS]

def diff_matrix(spectrum:Vector)->List[List[tuple]]: #returns a matrix of tuples: mass of possible neutral loss along with a computed tolerance
    return [ [(i - j, edge_tol(i,j)) if i > j else (0, 0) for i in spectrum] for j in spectrum]

def determine_losses(diff_matrix:List[List[tuple]])->List[List[Vector]]: #returns matrix of unambiguous neutral losses, as compomers
    dim = len(diff_matrix)
    mat = [[0]*dim for e in range(dim)]
    for i in range(dim):
        for j in range(dim):
            entry = diff_matrix[i][j]
            formulae = mf.HR_CF(alpha, entry[0], entry[1],use_heuristics=False, as_vec=True)
            if len(formulae) != 1:
                mat[i][j] = []
            else:
                mat[i][j] = formulae[0]
    return mat 


def edge_lst(spectrum:Vector)->List[List[tuple]]: #returns an edge list from the loss matrix
    loss_mat = determine_losses(diff_matrix(spectrum))
    return [(i,j) for i in range(len(loss_mat)) for j in range(len(loss_mat)) if loss_mat[i][j] != []]


######################################core algorithm######################################


def longest_path(spectrum:Vector)->tuple: #finds the longest path (to the molecular ion) through the subformula graph by chaining unambiguous losses
    res = []
    loss_mat = determine_losses(diff_matrix(spectrum))
    G = nx.DiGraph()
    G.add_edges_from(edge_lst(spectrum))
    p = nx.dag_longest_path(G)
    for i in range(len(p)-1):
        v1, v2 = p[i], p[i+1]
        res.append(loss_mat[v1][v2])
    return functools.reduce(mf.vec_add, res), p

#starting from a fragment ion, propose compomers within some tolerance. Then, these compomers are "fed forward" by adding the longest path compomer
#to the fragment ion compomers
def feed_forward(spectrum:Vector)->List[Vector]: 
    vec, path = longest_path(spectrum)
    frags = mf.HR_CF(alpha, spectrum[path[0]], vertex_tol(spectrum[path[0]]), as_vec=True)
    return [mf.vec_add(vec, i) for i in frags] 

#A list of candidates is proposed for each peak and each candidate formula for the molecular ion
#By constraining all other formulas to be subformulas of the chosen candidate
def feed_back(spectrum:Vector)->List[List[Vector]]: 
    res = []
    for f in feed_forward(spectrum):
        candidate_lst = [mf.HR_CF(alpha, i, vertex_tol(i,ppm=20),as_vec=True,parent=f) if mf.HR_CF(alpha, i, vertex_tol(i,ppm=20),as_vec=True,parent=f) != [] else [[]] for i in spectrum]
        res.append([(i, c[0]) for i, c in enumerate(candidate_lst)])
    return res


##########################subformula graph construction, scoring and visualisation##########################

def extract_f_list(candidates: List[List[tuple]])->List[List[Vector]]:
    return [c[1] for c in candidates if c[1] != []]


def extract_index(candidates: List[List[tuple]])->List[int]:
    return [c[0] for c in candidates if c[1] != []]


def subformula_graph(candidates:List[Vector], string=False)->nx.digraph: # generates the subformula graph for a given set of candidates
    G = nx.DiGraph()
    if not string:
         G.add_edges_from([ (i,j) for i in candidates for j in candidates if mf.is_subformula(i,j) and not mf.same_formula(i,j)])
    else:
        G.add_edges_from([(mf.produce_CF(i,alpha),mf.produce_CF(j,alpha)) for i in candidates for j in candidates if mf.is_subformula(i,j) and not mf.same_formula(i,j)])
    return G


def score_graph(G:nx.digraph)->float: #simple edge count, seems to work! Can modify...

    return nx.number_of_edges(G)/nx.number_of_nodes(G)


def analyse_MS(spectrum:Vector): #generate and score/rank subformula graphs
    candidates = feed_back(spectrum)
    graphs = [(subformula_graph(extract_f_list(l),string=True), extract_index(l)) for l in candidates]
    graphs.sort(key=lambda x: score_graph(x[0]), reverse=True)
    return graphs 


def draw_graph(graph:nx.digraph):
    nx.draw(graph,with_labels=True)
    plt.show()
    return




#takes the data from analyse_MS, along with a vector mass spectrum to produce a formula annotation
def return_formulae(spectrum:Vector, g:List[tuple],to_dict=False)->List[tuple]:
    formulae = list(nx.nodes(g[0]))
    formulae.sort(key=lambda x: mf.mw(mf.mf_to_compomer(x,alpha),alpha))
    index = g[1]
    if to_dict:
        return {int(round(i[1],0)): i[0] for i in zip(formulae, [spectrum[i] for i in index])}
    return list(zip(formulae, [spectrum[i] for i in index]))


def plot_annotated_mass_spectrum(spectrum, graph):
    max_mass = int(spectrum[-1][0]) + 50
    d = {int(round(p[0],0)): p[1] for p in spectrum}
    mz = [i for i in range(max_mass)]
    a = [d[i] if i in d else 0 for i in range(max_mass)]
    plt.bar(mz, a, color ='blue',
        width = 0.4)
    
    fs = return_formulae(extract_mass(spectrum),graph,to_dict=True)
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


#m = nc.mass_spectra(1203,nc.data,intensity_cutoff=0.5,base_norm=True)
#g = analyse_MS(extract_mass(m))

#plot_annotated_mass_spectrum(m,g[0])




#############################scoring (unused)##############################
def log_prob(val:float, mean:float, std=0.001)->float:
    z = abs(mean-val)/std
    res = math.erf(z/math.sqrt(2))
    return math.log(res) if res != 0.0  else -1000

def vertex_score(v:Vector, exp_mass:float)->float:  #score a vertex (formula) based off exp mass (log of gaussian)
    m = mf.mw(v,alpha)
    return log_prob(m,exp_mass,std=exp_mass/50000)
 
def edge_mass(e):
    return mf.mw(mf.diff(e[0],e[1]),alpha)

def edge_score(e:Vector, exp_mass:float)->float: #score an edge of the frag tree based off exp mass
    std_m1,std_m2 = mf.mw(e[0],alpha)/50000, mf.mw(e[1],alpha)/50000
    std_tot = math.sqrt(std_m1**2 + std_m2**2) 
    m_diff = edge_mass(e)
    return log_prob(m_diff,exp_mass,std=std_tot)

def edge_formulae(e):
    return [i[0] for i in e]

def masses(e):
    return e[0][1] -  e[1][1]

def sum_edge_scores(v): #given a vertex list, generate all edges and returns sum of edge scores
    e = subformula_graph(v)
    return sum([edge_score(edge_formulae(i),masses(i)) for i in e])

def sum_vertex_scores(v):
    return sum([vertex_score(i[0],i[1]) for i in v])

def total_score(v):
    return sum_edge_scores(v) + sum_vertex_scores(v)
