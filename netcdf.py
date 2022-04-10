import netCDF4 as nc
import numpy as np
import utils 
from typing import List

class NetCDF4ScanData:

    def __init__(self,path):
        data = nc.Dataset(path)
        self.no_of_peaks = data['point_count'][:]
        self.ion_current = data['total_intensity'][:]
        self.times = data['scan_acquisition_time'][:]
        self.masses = data['mass_values'][:]
        self.intensities = data['intensity_values'][:]
        self.maxscan = len(self.times)
        

class Chromatogram(NetCDF4ScanData):
    
    def get_TIC_peaks(self,intensity_cutoff:int)->List[int]: #get scan numbers of peaks in mass spectrum
        peak_list = utils.cluster([ind for ind, val in enumerate(self.ion_current) if val > intensity_cutoff])
        peak_numbers = []
        for peak in peak_list:
            d = [self.ion_current[i] for i in peak]
            peak_numbers.append(peak[d.index(max(d))])
        return peak_numbers
        
    def make_raw_MS(self,scan_number:int): #retrieves a raw mass spectrum
        masses, intensities = utils.split_lst(self.masses,self.no_of_peaks), utils.split_lst(self.intensities,self.no_of_peaks)
        return MassSpectrum(masses[scan_number-1],intensities[scan_number-1])
    
    def background(self, max_m:int)->np.ndarray: #stores average intensity of each bin into array (peaks >> avg intensity)
        background = np.zeros((max_m,))
        for ind, m in enumerate(self.masses):
            background[int(round(m,0))] += self.intensities[ind]
        return background / self.maxscan


class MassSpectrum:

    def __init__(self, masses:List[float], intensities:List[float]):
        self.masses = masses
        self.intensities = intensities
        self.total_intensity = sum(intensities)
    
    def get_MS(self): #returns a mass spectrum as a list of pairs
        return list(zip(self.masses,self.intensities))
    
    def conventional_norm(self): #normalise a mass spectrum according to convention (most intense peak = 100)
        conventional_norm_intensity = [i *100 /max(self.intensities) for i in self.intensities]
        return MassSpectrum(self.masses, conventional_norm_intensity)
    
    def de_noise(self, background:np.ndarray, cutoff_factor:float, intensity_cutoff:float): #de-noises a mass spectrum according to a cutoff tolerance for noise
        mass, intensity = [],[]
        for ind, m in enumerate(self.masses):
            point = self.intensities[ind]
            if point > background[int(round(m,0))] * cutoff_factor and point > intensity_cutoff * self.total_intensity:
                mass.append(m) 
                intensity.append(point)
        return MassSpectrum(mass,intensity)
    




#f = os.path.dirname(sys.argv[0]) + r'\WF12011.cdf'

#data = Chromatogram(f)
#ms = data.make_raw_MS(1204)
#print(ms.de_noise(data.background(500),5,0.005).get_MS())
    
    


    

