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
            formulae = []
            entry = diff_matrix[i][j]
            if entry[0] < 100: #added break statement
                formulae = mf.HR_CF(alpha, entry[0], entry[1],use_heuristics=False, as_vec=True)
            if len(formulae) != 1:
                mat[i][j] = []
            else:
                mat[i][j] = formulae[0]
    print("Finished computing loss matrix!")
    return mat 


def edge_lst(loss_mat:List[List[tuple]])->List[List[tuple]]: #returns an edge list from the loss matrix
    return [(i,j) for i in range(len(loss_mat)) for j in range(len(loss_mat)) if loss_mat[i][j] != []]


######################################core algorithm######################################


def longest_path(spectrum:Vector)->tuple: #finds the longest path (to the molecular ion) through the subformula graph by chaining unambiguous losses
    res = []
    loss_mat = determine_losses(diff_matrix(spectrum))
    G = nx.DiGraph()
    G.add_edges_from(edge_lst(loss_mat))
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

