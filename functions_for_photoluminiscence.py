# -*- coding: utf-8 -*-
"""
Auxiliary functions of "iterate_temperature_all_steps.py"
and for all processing steps

Mariano Barella

Dec 12th 2019

CIBION, Buenos Aires, Argentina
"""

import numpy as np
from scipy.optimize import curve_fit
import os
import json

## Uncomment these lines to ignore any warning message
import warnings
warnings.filterwarnings("ignore")

## DETERMINACION DE PARAMETROS
h = 4.135667516e-15 # in eV*s
c = 299792458 # in m/s

err_w0 = 0.005 # in um
instrumental_err = 0.05 # power meter uncertainty
skip_neighbours = 1 # number of spectrums to skip, for instance, first neighbour
skip_neighbours_threshold = 0.1 # relatvie value for excluding spectrum from the quotient

def closer(x,value):
    # returns the index of the closest element to value of the x array
    out = np.argmin(np.abs(x-value))
    return out

def manage_save_directory(path, new_folder_name):
    # Small function to create a new folder if not exist.
    new_folder_path = os.path.join(path, new_folder_name)
    if not os.path.exists(new_folder_path):
        os.makedirs(new_folder_path)
    return new_folder_path

def classification(value, totalbins, rango):
    # Bin the data. Classify a value into a bin.
    # totalbins = number of bins to divide rango (range)
    bin_max = totalbins - 1
    numbin = 0
    inf = rango[0]
    sup = rango[1]
    if value > sup:
        print('Value higher than max')
        return bin_max
    if value < inf:
        print('Value lower than min')
        return 0
    step = (sup - inf)/totalbins
    # tiene longitud totalbins + 1
    # pero en total son totalbins "cajitas" (bines)
    binned_range = np.arange(inf, sup + step, step)
    while numbin < bin_max:
        if (value >= binned_range[numbin] and value < binned_range[numbin+1]):
            break
        numbin += 1
        if numbin > bin_max:
            break
    return numbin

def lambda_to_energy(londa):
    # Energy in eV and wavelength in nm
    hc = 1239.84193 # Plank's constant times speed of light in eV*nm
    energy = hc/londa
    return energy

def energy_to_lambda(energy):
    # Energy in eV and wavelength in nm
    hc = 1239.84193 # Plank's constant times speed of light in eV*nm
    londa = hc/energy
    return londa

def lorentz2(x, *p):
    # Lorentz fitting function with an offset
    # gamma = FWHM
    # I = amplitude
    # x0 = center
    pi = 3.141592653589793
    I, gamma, x0, C = p
    return (1/pi) * I * (gamma/2)**2 / ((x - x0)**2 + (gamma/2)**2) + C

def temp_with_slope(Tzero, slope, Irrad):
    T = slope*Irrad + Tzero
    return T

def bose(energy, Temp):
    k = 0.000086173 # Boltzmann's constant in eV/K
    El = 2.3305 # excitation wavelength in eV
    aux = (energy - El) / (k*Temp)
    y = 1/( np.exp(aux) - 1 )
    return y

def log_bose(energy, cte, Temp):
    y = cte + np.log10(bose(energy, Temp))
    return y

def quotient_with_slope(energy, Tzero, slope, Irrad_n, Irrad_m):
    # Antistokes quotient fitting function
    # E is the energy array
    k = 0.000086173 # Boltzmann's constant in eV/K
    El = 2.3305 # excitation wavelength in eV
    Tn = temp_with_slope(Tzero, slope, Irrad_n)
    Tm = temp_with_slope(Tzero, slope, Irrad_m)
    aux_n = (energy - El) / ( k*Tn )
    aux_m = (energy - El) / ( k*Tm )
    # two equivalent expressions can be used for the quotient
    # reduced form:
#    quotient = A * np.exp( -(aux_n - aux_m) ) * np.sinh(aux_m) / np.sinh(aux_n)
    # explicit form:
    bose_m = 1/( np.exp(aux_m) - 1 )
    bose_n = 1/( np.exp(aux_n) - 1 )
    quotient = (Irrad_n * bose_n) / ( Irrad_m * bose_m )
    return quotient

def fit_lorentz2(p, x, y):
    return curve_fit(lorentz2, x, y, p0 = p)

def fit_bose(p, x, y):
    bose_func = lambda Temp : bose(x, Temp)
    return curve_fit(bose_func, x, y, p0 = p, bounds=(0, 5000))

def fit_log_bose(p, x, y):
    return curve_fit(log_bose, x, y, p0 = p)
    
def fit_quotient_for_beta(p, x, Tzero, Irrad_n, Irrad_m, y):
    quotient_lambda_func = lambda energy, beta : quotient_with_slope(energy, Tzero, beta, Irrad_n, Irrad_m)
    return curve_fit(quotient_lambda_func, x, y, p0 = p, bounds=(0, 5000))

def fit_quotient_for_Tzero(p, x, beta, Irrad_n, Irrad_m, y):
    quotient_lambda_func = lambda energy, Tzero : quotient_with_slope(energy, Tzero, beta, Irrad_n, Irrad_m)
    return curve_fit(quotient_lambda_func, x, y, p0 = p, bounds=(0, 5000))

def s2_single_quotient(init_params, energy, quotient_exp_data, Irrad_n, Irrad_m):
    # sum of squared residuals
    # function that is going to be minimized
    Tzero = init_params[0]
    beta = init_params[1]
    s2 = 0
    func = quotient_with_slope(energy, Tzero, beta, Irrad_n, Irrad_m)
    residuals = quotient_exp_data - func
    s2_residuals = residuals**2
    s2 = np.sum(s2_residuals)
    return s2

def calc_r2(observed, fitted):
    # Calculate coefficient of determination
    avg_y = observed.mean()
    # sum of squares of residuals
    ssres = ((observed - fitted)**2).sum()
    # total sum of squares
    sstot = ((observed - avg_y)**2).sum()
    return 1.0 - ssres/sstot

def calc_chi_squared(observed, expected):
    # Chi-squared (pearson) coefficient calculation
    aux = (observed - expected)**2/expected
    ans = aux.sum()
    return ans

#def gauss(x, A, mu, sigma):
#    return A*np.exp(-(x-mu)**2/(2*sigma**2))
#
#def fit_gauss(x, A, mu, sigma):
#    gauss_func = lambda energy, beta : quotient_with_slope(energy, Tzero, beta, Irrad_n, Irrad_m)
#    return curve_fit(gauss_func, x, y, p0 = p, bounds=(0, 5000))

def save_parameters(path_to, totalbins, meas_pow_bfp, Tzero, beta_guess, window, deg, repetitions ,\
                mode, factor, image_size_px, image_size_um, camera_px_length, \
                start_notch, end_notch, start_power, end_power, start_spr, max_spr, plot_flag, \
                lower_londa, upper_londa, alpha, R2th, single_NP_flag, \
                NP_int, do_step0, do_step1, do_step2, do_step3, do_step4):
    # save all input parameters of the analysis
    
    if not os.path.exists(path_to):
        os.makedirs(path_to)
        
    filepath = os.path.join(path_to,'input_parameters.txt')
    
    dictionary = {'totalbins':totalbins,
                  'meas_pow_bfp':meas_pow_bfp,
                  'Tzero':Tzero,
                  'beta_guess':beta_guess,
                  'window':window,
                  'deg':deg,
                  'repetitions':repetitions,
                  'mode':mode,
                  'factor':factor,
                  'image_size_px':image_size_px,
                  'image_size_um':image_size_um,
                  'camera_px_length':camera_px_length,
                  'start_notch':start_notch,
                  'end_notch':end_notch,
                  'start_power':start_power,
                  'end_power':end_power,
                  'start_spr':start_spr,
                  'max_spr':max_spr,
                  'plot_flag':plot_flag,
                  'lower_londa':lower_londa,
                  'upper_londa':upper_londa,
                  'alpha':alpha,
                  'R2th':R2th,
                  'single_NP_flag':single_NP_flag,
                  'NP_int':NP_int,
                  'do_step0':do_step0,
                  'do_step1':do_step1,
                  'do_step2':do_step2,
                  'do_step3':do_step3,
                  'do_step4':do_step4
                  }
    
    with open(filepath, 'w+') as f:
         f.write(json.dumps(dictionary))
    
    return

