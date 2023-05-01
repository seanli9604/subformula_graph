import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

from formula import Formula
from mass_spectrum import ExpMass, MassSpectrum
import spectrum_reader   


neutralised_casmi_spectrum = {124.01557645207001: 14.245388324299881, 126.01867645207: 0.7057755526892212, 138.01907645207: 4.123075377401849, 154.01407645207001: 14.469522232155256, 184.01207645207: 100.0}
example = MassSpectrum(neutralised_casmi_spectrum, 1)

harder_example = spectrum_reader.read_ms_from_file("MSBNK-CASMI_2016-SM878701.txt").read_mass_spectrum()
harder_example = MassSpectrum(harder_example, 1)

hardest_example = spectrum_reader.read_ms_from_file("MSBNK-CASMI_2016-SM878353.txt").read_mass_spectrum()
hardest_example = MassSpectrum(hardest_example, 1)

field_example = spectrum_reader.read_ms_from_file("example_chromatogram.CDF").get_baseline_subtracted_mass_spectrum(1203)
field_example = MassSpectrum(field_example, 10).conventional_norm()


def test_exp_mass():
    m1, m2 = ExpMass(205.113, 10), ExpMass(102.104, 20)
    sum_mass = m1 + m2 
    assert(sum_mass.mass == 307.217)
    formulae = [str(f) for f in m2.possible_formulae(["C","H", "N", "O"])]
    assert("C6H14O" == formulae[0])


def test_mass_spectrum():
    ms = MassSpectrum({23.1293123: 1, 54.1293192: 2, 59.12983: 10}, 0)
    new_ms = ms.conventional_norm()
    assert(new_ms.spectrum_dict == {23.1293123: 10, 54.1293192: 20, 59.12983: 100})
    new_ms = ms.get_p_filtered_mass_spectrum(54.1293192, 60)
    assert(new_ms.spectrum_dict == {54.1293192: 100})


def test_spectral_annotations():
    f = Formula.from_string("C6H4N2O5")
    frags = example.generate_fragment_formulae(f, delta_frag=1)
    assert(len(frags) == 3)
    frags = example.generate_fragment_formulae(f, delta_frag=4)
    assert(len(frags) == 5)
    annotations = example.get_spectral_annotations(["C", "H", "N", "O"], example.product_scoring_function, delta_frag=4)
    assert(len(annotations)) == 1
    t = annotations[0]
    assert(str(t) == "[C6H4O3, C5H4NO3, C6H4NO3, C6H4NO4, C6H4N2O5]")
    assert(example.vertex_scoring_function(t) == 5)
    assert(example.edge_scoring_function(t) == 9)
    assert(example.product_scoring_function(t) == 0.9)
    annotations = example.compute_most_likely_molecular_ion(["C","H","N", "O"], example.product_scoring_function, 125, 200, delta_frag=4)
    assert(str(annotations) == "[[C6H4O3, C5H4NO3, C6H4NO3, C6H4NO4, C6H4N2O5]]")
    
    annotations = harder_example.get_spectral_annotations(['C', 'H', 'N', 'O', 'S'], harder_example.product_scoring_function, delta_frag=1)
    assert(str(annotations[0][-1]) == "C13H10N2O3S" )
    assert([harder_example.edge_scoring_function(a) for a in annotations] == [47,6,3,1,1,0])
    assert([harder_example.product_scoring_function(a) for a in annotations] == [0.8545454545454545, 0.36363636363636365, 0.2727272727272727, 0.18181818181818182, 0.18181818181818182, 0])

    
    annotations = hardest_example.get_spectral_annotations(['C', 'H', 'N', 'O', 'F', 'S', 'Cl'], hardest_example.product_scoring_function, delta_frag=1)
    assert(str(annotations[0][-1]) == "C12H4N4OF6SCl2")
    assert(len(annotations) == 148)


def test_field_test():

    results = field_example.compute_most_likely_molecular_ion(['C', 'H', 'N', 'O', 'S'], field_example.product_scoring_function, 130, 250, delta_frag=3)
    assert(len(results) == 131)
    assert(str(results[1][-1]) == "C11H18N2O")
    assert(str(results[0][-1]) == "C11H20N2O")

    results.sort(key=field_example.edge_scoring_function, reverse=True)
    assert(str(results[1][-1]) == "C11H18N2O")

    results.sort(key=field_example.vertex_scoring_function, reverse=True)
    t = [(a[-1], field_example.vertex_scoring_function(a)) for a in results]
    assert(str(t[:15]) == "[(C11H20N2O, 65), (C11H18N2O, 63), (C11H16N2O, 61), (C9H17N5, 59), (C10H16N2O, 58), (C9H15N5, 57), (C6H16N4O3, 53), (C9H14N2O, 51), (C8H13N5, 49), (C6H19N5S, 49), (C8H11N5, 48), (C5H18N6O2, 47), (C8H19NO4, 46), (C11H18N2, 45), (C11H16N2, 43)]")
    assert(str(results[1][-1]) == "C11H18N2O")



