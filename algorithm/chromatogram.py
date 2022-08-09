import netCDF4 
import scipy.signal as sig
import matplotlib.pyplot as plt 
import utils 
from typing import List

class NetCDF4ScanData:

    def __init__(self,path):
        data = netCDF4.Dataset(path)
        self.no_of_peaks = data['point_count'][:]
        self.ion_current = data['total_intensity'][:]
        self.times = data['scan_acquisition_time'][:]
        self.masses = data['mass_values'][:]
        self.intensities = data['intensity_values'][:]
        self.maxscan = data.dimensions["scan_number"].size


class Chromatogram(NetCDF4ScanData):

    def avg_TIC(self):
        return sum(self.intensities)/len(self.intensities)
    
    def get_TIC_peak_list(self, intensity_factor):
        return sig.find_peaks(self.ion_current,height=self.avg_TIC()*intensity_factor)[0]
        
    def get_raw_MS(self,scan_number:int): #retrieves a raw mass spectrum
        masses, intensities = utils.split_lst(self.masses,self.no_of_peaks), utils.split_lst(self.intensities,self.no_of_peaks)
        spec_list = [MassSpectrum(masses[i],intensities[i]) for i in range(self.maxscan)]
        return spec_list[scan_number-1]
    
    def get_raw_MS_range(self, begin, end):
        masses, intensities = utils.split_lst(self.masses,self.no_of_peaks), utils.split_lst(self.intensities,self.no_of_peaks)
        spec_list = [MassSpectrum(masses[i],intensities[i]) for i in range(self.maxscan)]
        return spec_list[begin:end]
    
    def intensity_change(self, scan_no, peak_width, min_mass, max_mass):
        neighbouring_spec = self.get_raw_MS_range(scan_no - peak_width, scan_no + peak_width + 1)
        centre_spectrum = neighbouring_spec[peak_width]
        intensity_dict = {mass: [] for mass in centre_spectrum.get_MS(2) if mass > min_mass and mass < max_mass}
        for mass in centre_spectrum.get_MS(2):
            if mass > min_mass and mass < max_mass:
                for spec in neighbouring_spec:
                    ref_spec = spec.get_MS(2)
                    if mass in ref_spec:
                        intensity_dict[mass].append(ref_spec[mass])
                    else:
                        intensity_dict[mass].append(0)
        return intensity_dict


class MassSpectrum:

    def __init__(self, masses:List[float], intensities:List[float]):
        self.masses = masses
        self.intensities = intensities
        self.total_intensity = sum(intensities)
        self.bp_intensity = max(intensities)
        self.max_mass = int(max(self.masses)) + 1
    
    def get_MS(self, round_to_dp=False): #returns mass spectrum as a dictionary with masses as keys
        if round_to_dp:
            return {round(m, round_to_dp): i for m, i in zip(self.masses,self.intensities)}
        else:
            return {m: i for m, i in zip(self.masses,self.intensities)}
    
    def binned_MS(self): #if multiple masses in spectrum has same nominal mass, pick the one with maximum intensity
        mass_bin = {i: [] for i in range(self.max_mass + 1)}
        spectrum = self.get_MS()
        for m in spectrum:
            mass_bin[int(round(m,0))].append(m)
        new_masses = [max(mass_bin[i], key=lambda x: spectrum[x]) for i in mass_bin if mass_bin[i] != []]
        return MassSpectrum(new_masses, [spectrum[m] for m in new_masses]) 
    
    def intensity_cutoff(self, factor, max_mass): #removes peaks below a certain intensity or mass
        d = self.get_MS()
        MS_dict = {mass: d[mass] for mass in d if d[mass] > self.bp_intensity*factor or mass > max_mass}
        return MassSpectrum([m for m in MS_dict], [MS_dict[m] for m in MS_dict])

    def conventional_norm(self): #normalise a mass spectrum according to convention (most intense peak = 100)
        conventional_norm_intensity = [i / max(self.intensities) for i in self.intensities]
        return MassSpectrum(self.masses, conventional_norm_intensity)
    
    def find_candidate_masses(self, lower, upper):
        d = self.get_MS()
        trimmed_masses = [m for m in self.masses if m > lower and m < upper]
        clusters = utils.cluster_continuous(trimmed_masses)
        return [max(c,key=lambda x: d[x]) for c in clusters]

    def cluster_fragments(self, cluster_size=5):
        res = [[] for i in range (int(self.max_mass/cluster_size)+1)]
        for m, i in zip(self.masses, self.intensities):
            res[int(m/cluster_size)].append((m, i))
        ms = [max(l,key=lambda x: x[1]) for l in res if l != []]
        return MassSpectrum([i[0] for i in ms], [i[1] for i in ms])

        