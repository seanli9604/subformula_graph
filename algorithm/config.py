import json

#os.chdir(r"C:/Users/Darrin DeYoung/Desktop/abcde/subformula-graphs/")


class Parameters:

    def __init__(self, alphabet: list, heuristic: dict, mass_range: list, intensity_factor: float, peak_sensitivity: float, ppm_error: float):
        self.heuristic = heuristic
        self.mass_range = mass_range
        self.intensity_factor = intensity_factor
        self.peak_sensitivity = peak_sensitivity
        self.ppm_error = ppm_error
        self.alphabet = alphabet


def parse_config(config_path):
    with open(config_path, 'r') as config_fp:
        config = json.load(config_fp)
        return Parameters(
            [obj['element'] for obj in config["allowed_elements"]],
            {
                obj['element']: {
                    "min": obj["min"],
                    "max": obj["max"]
                }
                for obj in config["allowed_elements"]
                if obj.keys().__contains__("min") and obj.keys().__contains__("max")
            },
            [config["molecular_mass_range"]["min"],
                config["molecular_mass_range"]["max"]],
            config["intensity_factor"],
            config["chromatogram_peak_sensitivity"],
            config["ppm_error"]
        )