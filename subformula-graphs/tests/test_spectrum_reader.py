import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import spectrum_reader 

def test_casmi():
    data = spectrum_reader.read_ms_from_file("massbank_spectra_positive.txt")
    spectrum = data.read_mass_spectrum()
    assert(len(spectrum) == 228)
    spec_list = sorted([mass for mass in spectrum])
    assert(spec_list[-1] == 393.2072 - 1.00782503207 + 0.00054858)

    data = spectrum_reader.read_ms_from_file("massbank_spectra_negative.txt")
    spectrum = data.read_mass_spectrum()
    assert(len(spectrum) == 9)
    spec_list = sorted([mass for mass in spectrum])
    assert(spec_list[-1] == 183.0048 + 1.00782503207 - 0.00054858)



def test_recetox():
    data = spectrum_reader.read_ms_from_file("Recetox_example_spectra.txt")
    spectrum = data.read_mass_spectrum()
    assert(len(spectrum) == 19)
    spec_list = sorted([mass for mass in spectrum])
    assert(spec_list[-1] == 254.09985 + 0.00054858)

def test_netcdf():
    data = spectrum_reader.read_ms_from_file("example_chromatogram.CDF")
    spectrum = data.read_mass_spectrum(1203)
    spec_lst = sorted([mass for mass in spectrum], reverse=True, key=lambda m: spectrum[m])
    assert(int(spec_lst[0]) == 138)

def test_csv():
    data = spectrum_reader.read_ms_from_file("example.csv")
    spectrum = data.read_mass_spectrum()
    assert(len(spectrum) == 3)

