'''
Author: Sean Li
Date: May 1, 2023
This file contains the code necessary to read mass spectra, from a variety of input formats, and return a SpectrumReader object,
which can return the corresponding mass spectra a {mass: intensity} dictionary, in addition to metadata regarding the mass spectra 
(e.g ion mode, analyte formula, instument type). Baseline subtraction (only for the netCDF format) and corrections for electron/proton
mass are also performed at this stage, so we can deal with the masses of neutralised formulae. 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as 
published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. 
If not, see <https://www.gnu.org/licenses/>.
'''

from abc import ABC
import netCDF4  
import os

PROTON_MASS = 1.00782503207
ELECTRON_MASS = 0.00054858

def flatten(list_of_lists): # flattens a list of lists
    res = []
    for sublist in list_of_lists:
        for elem in sublist:
            res.append(elem)
    return res


#split a list into parts dictated by the elements of the parts array
# split_lst([1,2,3,4,5,6,7,8], [3,5]) -> [[1, 2, 3], [4, 5, 6, 7, 8]]
def split_lst(lst, parts): 
    res = []
    i,j = 0,0
    while j < len(parts):
        res.append(lst[i:i+parts[j]])
        i, j = i + parts[j], j + 1
    return res


def neutralise_formula_mass(ion_mode, mass): #neutralise a molecular mass so that the formula mass matches it
        
    if ion_mode == "POSITIVE":
        return mass + ELECTRON_MASS - PROTON_MASS
    elif ion_mode == "NEGATIVE":
        return mass - ELECTRON_MASS + PROTON_MASS
    else:
        raise Exception("Ion mode not found!")
 


class SpectrumReader(ABC):


    def __init__(self, filename):
        self.filename = filename

        try:
            open(filename) 
        except FileNotFoundError:
            raise Exception("File not found, please make sure you have correctly entered your file name!")
    
    
    def get_line(self, attr): #get a particular line in a text file
        return [line for line in self.get_raw_data() if line.startswith(attr)][0]
    

    def get_index(self, attr): #get the index of a particular line in a text file
        data = self.get_raw_data()
        return data.index([line for line in data if line.startswith(attr)][0])


    def extract_attribute(self, attr): #Given the start of a line, extract the attribute from the line
        return NotImplementedError


    def get_raw_data(self): #returns the needed data
        raise NotImplementedError


    def read_mass_spectrum(self):
        raise NotImplementedError
    

    def read_metadata(self):
        raise NotImplementedError



class JCampReader(SpectrumReader):


    def get_raw_data(self):
        return [line[:-1] for line in open(self.filename)]
    

    def get_mass_numbers(self): 
        start, end = self.get_index("##PEAK TABLE"), self.get_index("##END")
        return self.get_raw_data()[start+1: end]
    

    def read_mass_spectrum(self):
        mass_spectrum = {}
        numbers = flatten([line.split(" ") for line in self.get_mass_numbers()])
        for entry in numbers:
            pair = entry.split(",")
            mass, abundance = float(pair[0]), float(pair[1])
            mass_spectrum[mass] = abundance
        return mass_spectrum
    

    def extract_attribute(self, attr): 
        target_line = self.get_line(attr)
        return target_line.split("=")[1].replace(" ", "")
    

    def read_metadata(self):
        metadata = {"File Format": "JCAMP-DX",
                    "Compound Name": self.extract_attribute("##TITLE"),
                    "Molecular Formula": self.extract_attribute("##MOLFORM")
                    }
        return metadata



class MassBankReader(SpectrumReader):
    

    def get_raw_data(self):
        return [line[:-1] for line in open(self.filename)]


    def get_mass_numbers(self):
        start, end = self.get_index("PK$PEAK"), self.get_index("/")
        return self.get_raw_data()[start+1: end]
    
    def get_annotated_peaks(self):
        start, end = self.get_index("PK$ANNOTATION"), self.get_index("PK$NUM_PEAK")
        raw_data = self.get_raw_data()[start+1: end]
        annotations = [line.split(" ")[2:4] for line in raw_data]
        return {float(mass): formula[:-1] for mass, formula in annotations}
    

    def ion_mode(self):
        ion_mode = self.get_line("AC$MASS_SPECTROMETRY: ION_MODE")
        return ion_mode.split(" ")[-1]


    def read_mass_spectrum(self):
        mass_spectrum = {}
        ion_mode = self.ion_mode()
        numbers = [line.split(" ") for line in self.get_mass_numbers()]
        for entry in numbers:
            mass, abundance = neutralise_formula_mass(ion_mode, float(entry[2])), float(entry[3])
            mass_spectrum[mass] = abundance
        return mass_spectrum
    
    
    def extract_attribute(self, attr):
        target_line = self.get_line(attr)
        return target_line.split(":")[1].replace(" ", "")

    
    def read_metadata(self):
        metadata = {"File Format": "MassBank Data File",
                    "Compound Name": self.extract_attribute("CH$NAME"),
                    "SMILES": self.extract_attribute("CH$SMILES"),
                    "Ion Mode": self.extract_attribute("AC$MASS_SPECTROMETRY: ION_MODE"),
                    "Molecular Mass": float(self.extract_attribute("CH$EXACT_MASS:")),
                    "Molecular Formula": self.extract_attribute("CH$FORMULA:"),
                    "Instrument Type": self.extract_attribute("AC$INSTRUMENT_TYPE:")
                    }
        return metadata


class RecetoxReader(SpectrumReader):


    def get_raw_data(self):
        return [line[:-1] for line in open(self.filename)]


    def get_mass_numbers(self):
        start, end = self.get_index("Num Peaks"), len(self.get_raw_data()) -1
        return self.get_raw_data()[start+1: end]
    

    def ion_mode(self):
        ion_mode = self.get_line("AC$MASS_SPECTROMETRY: ION_MODE")
        return ion_mode.split(" ")[-1]


    def read_mass_spectrum(self):
        mass_spectrum = {}
        numbers = [line.split("\t") for line in self.get_mass_numbers()]
        for entry in numbers:
            mass, abundance = float(entry[0]) + ELECTRON_MASS  , float(entry[1])
            mass_spectrum[mass] = abundance
        return mass_spectrum
    
    
    def extract_attribute(self, attr):
        target_line = self.get_line(attr)
        return target_line.split(":")[1].replace(" ", "")

    
    def read_metadata(self):
        metadata = {"File Format": "Recetox Data File",
                    "Compound Name": self.extract_attribute("NAME"),
                    "SMILES": self.extract_attribute("SMILES"),
                    "Molecular Mass": float(self.extract_attribute("PRECURSORMZ")),
                    "Molecular Formula": self.extract_attribute("FORMULA"),
                    "Instrument Type": self.extract_attribute("INSTRUMENTTYPE")
                    }
        return metadata



class NetCDFReader(SpectrumReader):


    def __init__(self, filename):

        super().__init__(filename)

        self.data = netCDF4.Dataset(self.filename)
        self.no_of_peaks = self.data['point_count'][:]
        self.masses = self.data['mass_values'][:]
        self.intensities = self.data['intensity_values'][:]
        self.chromatogram = self.data["total_intensity"][:]
        self.time_values = self.data["scan_acquisition_time"][:] / 60  #in minutes!
        self.maxscan = self.data.dimensions["scan_number"].size
    

    def get_raw_data(self):
        return self.data


    def read_mass_spectrum(self, scan_number):
        masses, intensities = split_lst(self.masses,self.no_of_peaks), split_lst(self.intensities, self.no_of_peaks)
        return {mass: intensity for mass, intensity in zip(masses[scan_number-1], intensities[scan_number-1])}
    

    def get_baseline_subtracted_mass_spectrum(self, scan_number, baseline_shift=5):
        masses, intensities = split_lst(self.masses,self.no_of_peaks), split_lst(self.intensities, self.no_of_peaks)
        precise_mass_dict = {round(mass, 2): mass for mass in masses[scan_number -1]} #dictionary mapping rounded masses to masses
        spectrum  = {round(mass,2): intensity for mass, intensity in zip(masses[scan_number-1], intensities[scan_number-1])}
        baseline = {round(mass,2): intensity for mass, intensity in zip(masses[scan_number-1 - baseline_shift], intensities[scan_number-1 - baseline_shift])}

        for mass in baseline:
            if mass in spectrum:
                spectrum[mass] -= baseline[mass]

        return {precise_mass_dict[mass]: spectrum[mass] for mass in spectrum if spectrum[mass] > 0}


class CSVReader(SpectrumReader):

    '''Note: assume CSV has no headers and no metadata!'''

    def read_mass_spectrum(self):
        mass_spectrum = {}
        file = open(self.filename, "r")
        for line in file:
            mass, abundance = line.split(",")
            mass_spectrum[float(mass)] = float(abundance[:-1])
        return mass_spectrum




def read_ms_from_file(filename):

    ''' Given a suitable data file, parses it and converts it 
        into a reader object '''
        

    if filename.endswith(".jdx"):
        return JCampReader(filename)
    elif filename.split("/")[-1].startswith("Recetox"):
        return RecetoxReader(filename)
    elif filename.endswith(".txt"):
        return MassBankReader(filename)
    elif filename.endswith(".CDF") or filename.endswith(".cdf"):
        return NetCDFReader(filename)
    elif filename.endswith(".csv"):
        return CSVReader(filename)
    else:
        raise Exception("Error, cannot read this type of file!")


def read_ms_from_directory(dir):

    ''' Reads all mass spectra data files in a given directory '''
    
    spectra = []
    for index, filename in enumerate(os.listdir(dir)):
        path = dir + "/" + filename 
        spectra.append(read_ms_from_file(path))
    return spectra
