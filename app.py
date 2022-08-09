from flask import Flask, send_file, request
import main_fast
import matplotlib.pyplot as plt
from algorithm.plots import plot_chromatogram, plot_mass_spectrum
from algorithm.subformula_graph import find_molecular_ion

app = Flask(__name__)

@app.route("/<data>/chromatogram")
def chromatogram(data):
    data = main_fast.open_file(data)
    peak_sensitivity = int(request.args.get('peak_sensitivity', default= 50))
    plot_chromatogram(data, peak_sensitivity)
    return send_file('chromatogram.png', mimetype='image/png')

@app.route('/<data>/spectrum')
def spectrum(data):
    data = main_fast.open_file(data)
    scan_number = int(request.args.get('scan_number', default= 1))
    intensity_factor = float(request.args.get('intensity_factor', default= 0.001))
    ms = data.get_raw_MS(scan_number).binned_MS().conventional_norm().intensity_cutoff(intensity_factor, 999)
    plot_mass_spectrum(ms)
    return send_file('spectrum.png', mimetype='image/png')
    
# lmao copilot generated this
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