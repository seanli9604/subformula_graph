from flask import Flask, send_file, request
import main_fast
import matplotlib.pyplot as plt
from algorithm.plots import plot_chromatogram, plot_mass_spectrum
from algorithm.subformula_graph import find_molecular_ion

app = Flask(__name__)


def get_chromatagram(name):
    return main_fast.open_file(name)

def set_chromatagram(name, data):
    with open(name, 'w') as f:
        f.write(data)

def delete_chromatagram(name):
    try:
        os.remove(name)
    except OSError:
        pass

def mass_spectrum(name, scan_number, intensity_factor):
    data = get_chromatagram(name)
    ms = data.get_raw_MS(scan_number).binned_MS().conventional_norm().intensity_cutoff(intensity_factor, 999)
    return ms

@app.route("/<data>/chromatogram", methods=["GET"])
def chromatogram(data):
    data = main_fast.open_file(data)
    peak_sensitivity = int(request.args.get('peak_sensitivity', default=50))
    plot_chromatogram(data, peak_sensitivity)
    return send_file('chromatogram.png', mimetype='image/png')


@app.route('/<data>/spectrum', methods=["GET"])
def spectrum(data):
    scan_number = int(request.args.get('scan_number', default=1))
    intensity_factor = float(request.args.get(
        'intensity_factor', default=0.001))
    ms = mass_spectrum(data, scan_number, intensity_factor)
    plot_mass_spectrum(ms)
    return send_file('spectrum.png', mimetype='image/png')


@app.route('/<data>/formulae', methods=["GET"])
def formulae(data):
    peak_sensitivity = int(request.args.get('peak_sensitivity', default=50))
    alphabet = request.args.get('alphabet', default='ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    ppm_error = float(request.args.get('ppm_error', default=1))
    heuristic = request.args.get('heuristic', default='all')
    mass_range = [float(request.args.get('mass_range_min', default=0)), float(request.args.get('mass_range_max', default=1000))]
    scan_number = int(request.args.get('scan_number', default=1))
    intensity_factor = float(request.args.get('intensity_factor', default=0.001))
    ms = mass_spectrum(data, scan_number, intensity_factor)
    ions = find_molecular_ion(ms, *mass_range, alphabet, ppm_error, heuristic)
    return str(ions)


@app.route('/<data>', methods=["PUT", "POST", "GET", "DELETE"])
def data_routes(data):
    if request.method == "PUT":
        set_chromatogram(data, request.data)
        return "OK"
    elif request.method == "POST":
        set_chromatogram(data, request.data)
        return "OK"
    elif request.method == "GET":
        return main_fast.open_file(data)
    elif request.method == "DELETE":
        delete_chromatagram(data)
        return "OK"
    else:
        raise NotImplementedError()


@app.route('/<data>', methods=['GET'])
def index(data):
    return f'''
    <h1>Welcome to the Mass Spectrometry Analysis API</h1>
    <p>This API is used to retrieve mass spectrometry data from a netCDF4 file.</p>
    <p>To retrieve a chromatogram, use the following URL:</p>
    <p>/{data}/chromatogram?peak_sensitivity=50</p>
    <p>To retrieve a mass spectrum, use the following URL:</p>
    <p>/{data}/spectrum?scan_number=1&intensity_factor=0.001</p>
    '''
