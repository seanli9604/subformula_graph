import netCDF4 as nc
import os 
import sys
import formula as mf
import numpy as np
import itertools
from typing import List

f = os.path.dirname(sys.argv[0]) + r'\WF12011.cdf'
data = nc.Dataset(f)

def print_variables(data):
    for var in data.variables.values():
        print(var)

#############################reading raw data#############################

def datapoints(d)->np.ndarray:
    return d['point_count'][:]

def ion_current(d)->np.ndarray:
    return d['total_intensity'][:]

def times(d)->np.ndarray:
    return d['scan_acquisition_time'][:]


#############################utility functions#############################

def avg(l:List[float])->float:
    return sum(l)/len(l)

def cluster(l:List[float])->List[List[float]]: #cluster together contiguous scan numbers [1,2,5,6] -> [[1,2],[5,6]]
    res = [[l[0]]]
    for i in l:
        if i == res[-1][-1] + 1:
            res[-1].append(i)
        else:
            res.append([i])
    return res[1:]

#split a list into parts dictated by the elements of the parts array
# split_lst([1,2,3,4,5,6,7,8], [3,5]) -> [[1, 2, 3], [4, 5, 6, 7, 8]]
def split_lst(lst:List[float],parts:List[int])->List[List[float]]: 
    res = []
    i,j = 0,0
    while j < len(parts):
        res.append(lst[i:i+parts[j]])
        i, j = i + parts[j], j + 1
    return res


def re_norm(spectrum:List[tuple])->List[tuple]: #normalise spectrum according to convention (base peak = 100)
    max_int = max([i[1] for i in spectrum])
    f = lambda x: [x[0], x[1]*100/max_int]
    return [f(i) for i in spectrum]


def scans_from_time(time_range:tuple)->tuple: #retrieve scan range from time range
    l = [abs(time_range[0]*60-i) for i in times(data)]
    lower = l.index(min(l))
    u = [abs(time_range[1]*60-i) for i in times(data)]
    upper = u.index(min(u))
    return lower, upper

#############################chromatogram processing#############################

def locate_peaks(cutoff:float,chromatogram:List[float])->List[List[float]]: #locates and clusters peak given a cutoff definition of peak
    return cluster([i for i in range(len(chromatogram)) if chromatogram[i] > cutoff])

def separate_ms(data)->List[tuple]: #separate data into mass,intensity lists associated with each scan
    mass = split_lst(data['mass_values'][:], datapoints(data))
    intensity = split_lst(data["intensity_values"][:], datapoints(data))
    return [(mass[i],intensity[i]) for i in range(len(mass))]


def make_intensity_matrix(separated_data:List[tuple],lower_bound=45,upper_bound=500)->List[List[float]]: #actually make the intensity matrix
    res = np.zeros((upper_bound-lower_bound, len(separated_data)))
    for i in range(len(separated_data)):  #M_xy, where x corresponds to mass and y corresponds to scan number. Entry corresponds to intensity
        for j in range(len(separated_data[i][0])):
            ind = round(separated_data[i][0][j]) - lower_bound
            res[ind,i] += separated_data[i][1][j]
    return res  


#returns mass spectrum for a given scan number, rudimentary de-noising
def mass_spectra(scan_no:int, data, intensity_cutoff=0, base_norm=False)->List[tuple]: 
    mat = make_intensity_matrix(separate_ms(data))
    norm = ion_current(data)[scan_no-1]
    p = (sum(data['point_count'][:scan_no-1]), sum(data['point_count'][:scan_no]))
    mass = data['mass_values'][p[0]:p[1]]
    raw_intensity = data['intensity_values'][p[0]:p[1]]
    normed_intensity = [i*100/norm for i in raw_intensity]
    
    res = []
    for i in range(len(mass)):
        bin = int(round(mass[i],0)) - 45
        if raw_intensity[i] > avg(mat[bin,:]) * 2 and normed_intensity[i] > intensity_cutoff:
            res.append([mass[i],normed_intensity[i]])

    return res if not base_norm else re_norm(res)

##################################retrieving ms#################################


def find_MS(data, ic_cutoff:int)->List[int]: #retrieve peak index associated with maximum-intensity scan, rather than average
    peak_lst = locate_peaks(ic_cutoff, ion_current(data))
    res = []
    for peak in peak_lst:
        d = [ion_current(data)[i] for i in peak]
        res.append(peak[d.index(max(d))])
    return res


def bin_masses(peak, data): #given list of spectra, bin masses together and average them
    spec_lst = [ [j[0] for j in mass_spectra(i,data)] for i in peak]
    masses = list(itertools.chain(*spec_lst))
    masses.sort()
    return [avg(list(g)) for k,g in itertools.groupby(masses,key=lambda x: round(x,0))]

