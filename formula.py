import itertools
import functools
import re
import numpy as np
from typing import List, TypeVar


tree = TypeVar('Tree')

#these are all for the most abundant isotope!
int_masses = {"H": 1,
                  "C": 12,
                  "N": 14,
                  "O": 16,
                  "F": 19,
                  "P": 31,
                  "S": 32,
                  "Cl": 35}

exact_masses =  {"H": 1.00782503207,
                  "C": 12,
                  "N": 14.0030740048,
                  "O": 15.99491461956,
                  "F": 18.99840322,
                  "P": 30.97376163,
                  "S": 31.97207100,
                  "Cl": 34.96885268}

mass_of_electron = 0

class Formula:
    
    def __init__(self, alphabet:List[str], compomer:List[int]):
        self.alphabet = alphabet
        self.compomer = np.array(compomer)
        self.nom_mass = sum(map(lambda x, y: int_masses[x]*y, self.alphabet, self.compomer))
        self.exact_mass = sum(map(lambda x, y: exact_masses[x]*y, self.alphabet, self.compomer))
    
    def __repr__(self):
        
        def func(elem, coeff):
            if coeff == 0:
                return ""
            elif coeff == 1:
                return elem
            else:
                return elem + str(coeff)

        partial_str = map(func, self.alphabet, self.compomer)
        return functools.reduce(lambda x,y: x + y, partial_str)
    

    def formula_from_string(string:str, alphabet: List[str]): #could rewrite this
        compomer = [0]*len(alphabet)
        for num, elem in enumerate(alphabet):
            if re.search(elem + "[0-9][0-9]",string):
                ind = re.search(elem + "[0-9][0-9]",string).end()
                compomer[num] += int(string[ind-2:ind])
            elif re.search(elem + "[0-9]",string):
                ind = re.search(elem + "[0-9]",string).end()
                compomer[num] += int(string[ind-1])
            elif re.search(elem,string):
                compomer[num] += 1
        return Formula(alphabet,compomer)
    
    def add_formulas(f1, f2):
        return Formula(f1.alphabet, f1.compomer + f2.compomer)
    
    def is_subformula(formula, sub):
        return min(formula.compomer - sub.compomer) >= 0 
    
    def same_formula(f1,f2):
        return np.linalg.norm(f1.compomer - f2.compomer) == 0
    
    def subtract_formulas(f1, f2):
        return Formula(f1.alphabet, f1.compomer - f2.compomer)
    
    def mass_deviation(self, exp_mass): #return mass deviation from experimental mass in ppm
        return 10**6 * (self.exact_mass - exp_mass) / exp_mass
    
    def dbe(self): #ring + double bond formula: DBE = n(C) + 1 - n(H)/2 + n(N)/2 
        return 1 + self.compomer[self.alphabet.index("C")] - (self.compomer[self.alphabet.index("H")]/2) + (self.compomer[self.alphabet.index("N")]/2)


class Tree:
    
    def __init__(self, data=None, daughters=None):
        self.daughters = daughters
        self.data = data

    def traverse(self, index):
        return self.daughters[index]
    
    def sum_tree(self):
        if self.data == None:
            return 0
        elif self.daughters == None:
            return self.data 
        else:
            return self.data + self.traverse(0).sum_tree()
    

    def expand_tree(self, alphabet): #expands formula tree into list of compomers with an alphabet
        res = []

        def s(t, r , depth):
            if depth > 0:
                r[depth-1] = t.data 
            if depth == len(alphabet):
                res.append(r.copy())
            else:
                for daughter in t.daughters:
                    s(daughter, r, depth+1)

        s(self, np.zeros((len(alphabet),),dtype=int), 0)
        return [i[::-1] // alphabet for i in res]


class ElemRange:

    def __init__(self, alpha, mass):
        self.alpha = alpha
        self.mass = mass
        self.min_elems = {elem: 0 for elem in alpha}
        self.max_elems = {elem: int(mass / int_masses[elem]) for elem in alpha}
    
    def set_min_max_elements(self, heuristic_list): 
        for elem, entry in heuristic_list.items():
            if "min" in entry:
                self.min_elems[elem] = entry['min']
            if "max" in entry:
                self.max_elems[elem] = entry['max'] 
        return

    def default_heuristic(self):
        heuristic = {"C": {"min": int(self.mass/(4*12))}, "H": {"max": 2 * self.max_elems["C"] + 2}}
        self.set_min_max_elements(heuristic)
        return
    
    def custom_heuristic(self, heuristic):
        res = {"C": {"min": int(self.mass/(4*12))}, "H": {"max": 2 * self.max_elems["C"] + 2}}
        for elem in heuristic:
            res[elem] = heuristic[elem]
        self.set_min_max_elements(heuristic)
        return

    def parental_restriction(self, formula): #fragments must be subformulae of parent!
        f = {elem : coeff for elem,coeff in zip(formula.alphabet,formula.compomer)}
        for elem, coeff in f.items():
            self.max_elems[elem] = coeff
        return

    def get_gfs(self):
        gfs = []
        for elem in self.alpha:
            gf = range(int_masses[elem] * self.min_elems[elem], int_masses[elem] * self.max_elems[elem] + 1, int_masses[elem])
            gfs.append(list(gf))
        return gfs


class MFGenerator(ElemRange):

    def enumerate_compomers(self):

        def build_trees_from_treelist(generating_fn:List[int], list_of_tree_lists:List[List[tree]])->List[tree]:
            first_tree = lambda p: p[1][0]
            return [Tree(prod[0], prod[1]) for prod in itertools.product(generating_fn, list_of_tree_lists) 
            if first_tree(prod).sum_tree() + prod[0] < self.mass + 1]
    
        def bin_trees(tree_list:List[tree])->List[List[tree]]:
            res = [[] for i in range(self.mass + 1)]
            for tree in tree_list:
                if tree.sum_tree() < self.mass + 1:
                    res[tree.sum_tree()].append(tree)
            return [tuple(trees) for trees in res if len(trees) != 0]

        base_case = [[Tree()]] #singleton list of lists
        list_of_tree_lists = functools.reduce(lambda t,gf: bin_trees(build_trees_from_treelist(gf,t)), self.get_gfs(), base_case)
        
        formula_tree = Tree(0, list_of_tree_lists[-1]) #merge final list of trees into one tree with dummy node
        alpha_array = np.array([int_masses[elem] for elem in self.alpha],dtype=int)

        return formula_tree.expand_tree(alpha_array)
    
    
    def get_formula_list(self, exp_mass, tol, DBE_restriction = None): #e.g lambda x: x.dbe() >= 0 and x.dbe().is_integer()
        formula_list = [Formula(self.alpha, c) for c in  self.enumerate_compomers()]
        formula_list = filter(lambda x: abs(x.mass_deviation(exp_mass)) < tol, formula_list)
        if DBE_restriction != None:
            formula_list = filter(DBE_restriction, formula_list)
        return sorted(list(formula_list), key=lambda x: x.mass_deviation(exp_mass))
