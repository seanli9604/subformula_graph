import algorithm.chromatogram as netcdf
import matplotlib.pyplot as plt 

def plot_chromatogram(data:netcdf.Chromatogram, intensity_factor:int):
    y = data.ion_current
    t = [data.times[i]/60 for i in range(data.maxscan)]
    plt.plot(t,y)

    peaks = data.get_TIC_peak_list(intensity_factor)

    for p in peaks:
        plt.annotate("{} \nRT={}".format(p, round(t[p],2)), # this is the text
                (t[p], y[p]), # these are the coordinates to position the label
                fontsize=8,
                textcoords="offset points", # how to position the text
                xytext=(0,2), # distance from text to points (x,y)
                ha='center') # horizontal alignment can be left, right or center

    plt.ylabel("Total ion current")
    plt.xlabel("Time (minutes)")
    plt.show()


def plot_mass_spectrum(spectrum:netcdf.MassSpectrum):
    max_mass = int(max(spectrum.masses)) + 1
    ms_info = spectrum.get_MS()
    ms_info_rounded = {int(round(m,0)): ms_info[m] for m in ms_info}


    mz = [i for i in range(max_mass)]
    a = [ms_info_rounded[i] if i in ms_info_rounded else 0 for i in range(max_mass)]

    plt.bar(mz, a, color ='blue',
        width = 0.4)
    
    for m in spectrum.find_candidate_masses(50, max_mass):
        plt.annotate(round(m,4)  , # this is the text
                (int(round(m,0)), ms_info[m]),
                fontsize=8, 
                textcoords="offset points", # how to position the text
                xytext=(0,2), # distance from text to points (x,y)
                ha='center') # horizontal alignment can be left, right or center

    plt.xlabel("m/z")
    plt.ylabel("Abundance")
    plt.title("Mass spectrum")
    plt.show()


def plot_annotated_mass_spectrum(spectrum, graph): 
    max_mass = int(spectrum[-1][0]) + 50
    d = {int(round(p[0],0)): p[1] for p in spectrum}
    mz = [i for i in range(max_mass)]
    a = [d[i] if i in d else 0 for i in range(max_mass)]
    plt.bar(mz, a, color ='blue',
        width = 0.4)
    
    #fs = ft.return_formulae(ft.extract_mass(spectrum),graph,to_dict=True)
    for p in spectrum:
        if int(round(p[0],0)) in fs:
            plt.annotate(fs[int(round(p[0],0))], # this is the text
                (int(round(p[0],0)), p[1]),
                fontsize=8, 
                textcoords="offset points", # how to position the text
                xytext=(0,2), # distance from text to points (x,y)
                ha='center') # horizontal alignment can be left, right or center

    plt.xlabel("m/z")
    plt.ylabel("Abundance")
    plt.title("Mass spectrum")
    plt.show()
