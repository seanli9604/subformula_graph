import os
import re 

#os.chdir(r"C:/Users/Darrin DeYoung/Desktop/abcde/subformula-graphs/")

class Parameters:

    def __init__(self, alphabet:list, heuristic:dict, mass_range:list, intensity_factor:float, peak_sensitivity:float, ppm_error:float):
        self.heuristic = heuristic 
        self.mass_range = mass_range
        self.intensity_factor = intensity_factor
        self.peak_sensitivity = peak_sensitivity
        self.ppm_error = ppm_error
        self.alphabet = alphabet

def parse_line(line, pattern, func): #parses a line according to specified regex pattern and function which typecasts the match
    patt = re.compile(pattern) 
    match = patt.search(line)[0]
    return func(match)

def make_heuristic_list(tokens):
    alphabet = []
    heuristic_list = {}
    patt = re.compile("[0-9]+")
    p2 = re.compile("[CHNOPSF(Cl)(Si)]")
    for t in tokens:
        elem = p2.search(t)[0]
        alphabet.append(elem)
        nums = patt.findall(t)
        if nums != []:
            heuristic_list[elem] = {"min": int(nums[0]), "max": int(nums[1])}
    return (alphabet, heuristic_list)


def parse_main():
    f = "".join([line for line in open("config.txt")])
    f = f.split("}")
    parameters = {}

    for l in f:

        if "allowed elements" in l.lower():
            l = l.split(":")
            f = lambda s: make_heuristic_list(s.replace(" ", "").split(","))
            parameters["alphabet"], parameters["elem"] = parse_line(l[1], "[CHNOPSF(Cl)(Si)].*[CHNOPSF(Cl)(Si)]", f)
        
        elif "molecular mass range" in l.lower():
            f = lambda s: [int(i) for i in s.replace(" ", "").split(",")]
            parameters["mass_range"] = parse_line(l, "[0-9]+ *, *[0-9]+", f)
        
        elif "intensity factor" in l.lower():
            parameters["int_factor"] = parse_line(l, "0.[0-9]*", float)
        
        elif "chromatogram peak sensitivity" in l.lower():
            parameters["peak_sensitivity"] = parse_line(l, "[0-9]+", float)
        
        elif "ppm error" in l.lower():
            parameters["ppm err"] = parse_line(l, "[0-9]+", int)

    return Parameters(parameters["alphabet"], parameters["elem"],parameters["mass_range"], 
    parameters["int_factor"], parameters["peak_sensitivity"], parameters["ppm err"]) 

