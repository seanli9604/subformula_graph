import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

from formula import Formula
from mass_spectrum import ExpMass, MassSpectrum
import spectrum_reader
from visualisation import MSVisualiser  



harder_example = spectrum_reader.read_ms_from_file("MSBNK-CASMI_2016-SM878701.txt").read_mass_spectrum()
harder_example = MassSpectrum(harder_example, 1).conventional_norm()

field_example = spectrum_reader.read_ms_from_file("example_chromatogram.CDF").get_baseline_subtracted_mass_spectrum(1203)
field_example = MassSpectrum(field_example, 10).conventional_norm()


def test_vis():
    annotations = harder_example.get_spectral_annotations(['C', 'H', 'N', 'O', 'S'], harder_example.product_scoring_function, delta_frag=1)
    correct_annotation = annotations[0]
    vis = MSVisualiser(harder_example, correct_annotation)
    assert(vis.formula_dict == {'C5H5N': 0.15396195263205625, 'C6H5N': 0.17080983216662443, 'C6H5NO': 0.6864252915773101, 'C12H9N2': 0.32760251172884614, 'C13H8N2': 0.7663026847804278, 'C13H9N2': 13.77437400341001, 'C13H10N2': 3.6260218455459468, 'C13H9N2O': 0.43266530686329946, 'C13H10N2O': 2.2544509410181997, 'C13H9N2O2': 3.3937708573585312, 'C13H10N2O3S': 100.0})
    assert(str(vis.generate_formula_labels()) == '[C5H5N, C6H5N, C6H5NO, C12H9N2, C13H9N2, C13H10N2O, C13H9N2O2, C13H10N2O3S]')

    vis.plot_ms_as_2D_fragment()

    #results = field_example.compute_most_likely_molecular_ion(['C', 'H', 'N', 'O', 'S'], field_example.product_scoring_function, 130, 250, delta_frag=3)
    #correct_annotation = results[1]
    #vis = MSVisualiser(field_example, correct_annotation)
    #vis.plot_ms_as_2D_fragment()

    