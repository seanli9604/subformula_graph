import itertools
import re
from typing import List
from typing import Tuple

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

alphabet = [1,12,14,16,19,31,32,35]
exact_m = [1.00782503207,12,14.0030740048,15.99491461956,18.99840322,30.97376163,31.97207100,34.96885268]

#######################types###########################
Vector = List[float]


#######################formula/compomer manipulation###########################

def mf_to_compomer(mf:str, elements:List[str])->List[int]: #converts a string molecular formula into a compomer given a list of alphabet
    
    if mf == []:  #no parent
        return []

    compomer = [0]*len(elements)
    for num, elem in enumerate(elements):
        if re.search(elem + "[0-9][0-9]",mf):
            ind = re.search(elem + "[0-9][0-9]",mf).end()
            compomer[num] += int(mf[ind-2:ind])
        elif re.search(elem + "[0-9]",mf):
            ind = re.search(elem + "[0-9]",mf).end()
            compomer[num] += int(mf[ind-1])
        elif re.search(elem,mf):
            compomer[num] += 1
    return compomer

def alpha_to_vec(alpha: List[str],exact=False)->Vector:
    
    if exact:
        vec = [exact_masses[i] for i in alpha]
    else:
        vec = [int_masses[i] for i in alpha]
    return vec

def mw(compomer:Vector, alpha:List[str])->float:
    
    def dot_product(a:Vector,b:Vector)->float:
        return sum([a[i]*b[i] for i in range(len(a))])
    
    v = alpha_to_vec(alpha, exact=True)
    return dot_product(v, compomer) 

def vec_add(c1:Vector, c2:Vector)->Vector:
    return [c1[i] + c2[i] for i in range(len(c1))]

def is_subformula(f:Vector, s:Vector)->bool:
    return sum([f[i] - s[i] < 0 for i in range(len(f))]) == 0

def same_formula(f1:Vector,f2:Vector)->bool: 
    return sum([abs(f1[i]-f2[i]) for i in range(len(f1))]) == 0

def diff(f:Vector, s:Vector)->Vector: # Vector subtraction
    return [f[i] - s[i] for i in range(len(f))]

def mass_diff(formula, mass:float, alphabet:List[str], as_string=False): #mass deviation between "exact" mass of possible formula and measured mass
    return abs(mass  - mw(formula, alphabet)) if not as_string else abs(mass - mw(mf_to_compomer(formula,alphabet), alphabet))


#########################formula enumeration algorithm#########################


def enum_mf(alphabet: Vector, M:int, use_heuristics=True)->List[List[int]]:
    
    def is_leaf(t:tuple)->bool:
        return len(t) == 1 and type(t[0]) is int

    def sum_tree(t:tuple)->int:
        if is_leaf(t):
            return sum(t)
        elif type(t[0]) is int:
            return t[0] + sum_tree(t[1])
        else:
            return sum_tree(t[0])

    def bin_into(l:tuple,maxint:int)->List[tuple]: #bin subtrees into groups
        res = [[] for i in range(maxint)]
        for t in l:
            if sum_tree(t) < maxint:
                res[sum_tree(t)].append(t)
        return [tuple(i) for i in res if len(i) != 0]

    def traverse_tree(tree:tuple)->List[List[int]]: #traverse formula tree to obtain formula list
        res = []
        def s(t,r):
            if is_leaf(t):
                res.append(r[1:]+[t[0]])
            elif type(t[0]) is int:
                [s(branch,r+[t[0]]) for branch in t[1:]] 
            else:
                [s(branch, r) for branch in t]
        s(tree,[])
        return res

# Produce generating functions based on the alphabet
    def produce_gfs(alphabet:Vector, M:int, use_heuristics:bool)->List[List[int]]:
        l = [ [ j for j in range(0,M+1,i)] for i in alphabet]
        if use_heuristics:
            l[0] = l[0][:2* len(l[1]) + 2] # truncate hydrogens
            l[1] = l[1][int(M/(4*12)):] #truncate carbons
        return l 

    def tree(alpha:Vector, M:int)->tuple: #construct formula tree
        L = produce_gfs(alpha,M, use_heuristics)
        def r(L,M,prod):
            return list(prod)[-1] if L == [] else r(L[1:], M, bin_into(itertools.product(L[0],prod),M+1))
        return (0, r(L[1:],M,[(i,) for i in L[0]]))
    
    def compomer(vec:Vector)->Vector: #obtain c_i of compomer by n_i/m_i
        return [int(vec[i]/alphabet[i]) for i in range(len(vec))] 
        
    return [compomer(i[::-1]) for i in traverse_tree(tree(alphabet, M))]


#########################HR-MS functions#########################

#ad-hoc dbe function!
def dbe(compomer:Vector)->float:
    return compomer[0] + 1 - (compomer[1]/2) + (compomer[2]/2)

def HR_CF(alphabet:List[str], M:float, tol:float, parent=[], use_heuristics = True, as_vec = False)->List[List[int]]:
    lst = enum_mf(alpha_to_vec(alphabet), int(round(M,0)), use_heuristics)
    #filters list of formulae by accuracy of mass spectral data
    lst = [i for i in lst if mass_diff(i, M, alphabet) < tol ]
    lst.sort(key=lambda x: mass_diff(x, M, alphabet))
    if parent != []:
        lst = [i for i in lst if is_subformula(parent, i) and dbe(i) >= 0] 
    return lst if as_vec else CF(lst, alphabet)


def produce_CF(compomer:Vector, alphabet:List[str])->str: 
    # Can be any list of symbols corresponding to elements
    a = alphabet.copy()
    #a[0],a[1] = a[1],a[0]
    res = [a[i]+str(compomer[i]) for i in range(len(compomer)) if compomer[i] != 0]
    return  "".join([i if len(i) != 2 or i[1] != "1" else i[0] for i in res])

def CF(compomer_lst:Vector, alphabet:List[str])->List[str]: #returns list of molecular formulae
        return [produce_CF(i, alphabet) for i in compomer_lst]
