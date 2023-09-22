import os
import random 
import sys
sys.path.append('subformula-graphs')

import matplotlib.pyplot as plt

import spectrum_reader
from formula import Formula, FormulaGenerator
from mass_spectrum import MassSpectrum 
import constants
from visualisation import MSVisualiser

spectra_dir = "mass-spectra/CASMI_2016"

def ppm(mass1, mass2):
    return abs(mass1 - mass2) * 10**6 / mass1

def mix_spectra(spec1, spec2, ppm_criterion):
    '''Mix two mass spectra together into a combined mass spectra.
    (ppm_criterion dictates how close masses have to be before being combined)'''
    norm1, norm2 = sum(spec1.values()), sum(spec2.values())
    merged = {mass: intensity / norm1 for mass, intensity in spec1.items()}
    for m2, intensity in spec2.items():
        m1 = min(spec1.keys(), key=lambda m1: ppm(m1, m2))
        if ppm(m1, m2) < ppm_criterion:
            merged[m1] += intensity / norm2
        else:
            merged[m2] = intensity / norm2
    pairs = sorted(merged.items())
    return {mass: intensity for mass, intensity in pairs}


def sample_spectra(seed, num_samples, mass_threshold):
    '''Randomly sample (with replacement) pairs of mass spectra (containing the molecular ion), 
    given a seed value for the random number generator'''
    res = []
    spectra = os.listdir(spectra_dir)
    random.seed(seed)
    while len(res) < num_samples:
        pair = random.choices(spectra, k=2)
        if is_valid_spectra(pair[0], mass_threshold) and is_valid_spectra(pair[1], mass_threshold):
            res.append(pair)
    return res 


def is_valid_spectra(file, mass_threshold):
    '''Check if the spectra contained in a file is "valid" (i.e contains molecular ion, below a threshold mass
    and has more than one mass peak'''
    ms_file = spectrum_reader.read_ms_from_file(spectra_dir + "/" + file)
    spec = ms_file.read_mass_spectrum()
    formula = Formula.from_string(ms_file.read_metadata()["Molecular Formula"])
    closest_mass = min(spec.keys(), key=lambda m: abs(formula.exact_mass - m))
    within_mass_threshold = max(spec.keys()) < mass_threshold
    return abs(closest_mass - formula.exact_mass) < 0.1 and len(spec) > 1 and within_mass_threshold


def extract_relevant_info(pair, ppm_criterion):
    '''Get the combined mass spectra + metadata needed for testing given two file names'''
    file1 = spectrum_reader.read_ms_from_file(spectra_dir + "/" + pair[0])
    file2 = spectrum_reader.read_ms_from_file(spectra_dir + "/" + pair[1])
    f1 = Formula.from_string(file1.read_metadata()["Molecular Formula"])
    f2 = Formula.from_string(file2.read_metadata()["Molecular Formula"])
    defaultalpha = set(["C", "H", "N", "O"])
    defaultalpha = defaultalpha.union(set(f1.get_constituent_elements()))
    defaultalpha = list(defaultalpha.union(set(f2.get_constituent_elements())))
    defaultalpha.sort(key=lambda x: constants.EXACT_MASSES.get(x))

    mixed_spec = mix_spectra(file1.read_mass_spectrum(), file2.read_mass_spectrum(), ppm_criterion)
    ms = MassSpectrum(mixed_spec, 1)
    lower, upper = min(ms.masses), max(ms.masses) + 1
    print(f1.exact_mass, f2.exact_mass, upper)
    return f1, f2, defaultalpha, lower, upper, ms 


def separation_is_successful(cand, f1, f2):
    '''Checks if both cands are correctly identified'''
    dominant_compound = cand[0][-1]
    second_compound = None
    for c in cand:
        if not c[-1].is_subformula(dominant_compound):
            second_compound = c[-1]
            break
    
    if not second_compound:
        return False #exit early
    
    case1 = (f1 == dominant_compound and f2 == second_compound)
    case2 = (f2 == dominant_compound and f1 == second_compound)
    return (case1 or case2)


def test_combined_ms(seed, num_samples, mass_threshold, ppm_criterion = 10):
    '''Run test on mixed spectra'''
    res = []
    pairs = sample_spectra(seed, num_samples, mass_threshold)
    for p in pairs:
        f1, f2, defaultalpha, lower, upper, ms = extract_relevant_info(p, ppm_criterion)
        print(f1, f2)
        cand = ms.compute_most_likely_molecular_ion(defaultalpha, ms.product_scoring_function, lower, upper)
        r = separation_is_successful(cand, f1, f2)
        print(r)
        subformulae = f1.is_subformula(f2) or f2.is_subformula(f1)
        res.append((f1, f2, subformulae, r))
    
    subformulae_excluded_total = sum(not r[2] for r in res)
    total_separated = sum(r[3] for r in res)
    print(total_separated / num_samples)
    print(total_separated / subformulae_excluded_total)
    return res 

# seed for random number generator used: 999
# masses tested: 150, 200, 250, 300, 350, 400

# test_combined_ms(999, 100, 150)
