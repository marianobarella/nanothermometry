# -*- coding: utf-8 -*-
"""
Process of single AuNPs spectral confocal image 
acquired with PySpectrum at CIBION

Mariano Barella

16 aug 2019

based on "witec_data_photoluminiscense.py"

"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from functions_for_photoluminiscence import manage_save_directory, closer, \
    instrumental_err, err_w0

try:
    plt.style.use('for_confocal.mplstyle')
except:
    print('Pre-defined matplotlib style was not loaded.')

plt.ioff()
plt.close('all')

def process_irrad_calibration(folder, path_to, calibration_image_size_px, 
                               calibration_image_size_um, camera_px_length, window, 
                               deg, repetitions, mode, factor, meas_pow_bfp, 
                               start_notch, end_notch, start_power, end_power, start_spr, 
                               lower_londa, upper_londa, exp_time_fast, w0):
    # plot limits
    y_min_spec_cal = 8000
    y_max_spec_cal = 238000
    
    NP = folder.split('Spectrum_')[-1]
    
    save_folder = os.path.join(path_to, NP)
    
    common_path = os.path.join(path_to,'common_plots')
    
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
            
    list_of_files = os.listdir(folder)
    wavelength_filename = [f for f in list_of_files if re.search('wavelength',f)]
    list_of_files.sort()
    list_of_files = [f for f in list_of_files if not os.path.isdir(folder+f)]
    list_of_files = [f for f in list_of_files if ((not os.path.isdir(folder+f)) \
                                                  and (re.search('_i\d\d\d\d_j\d\d\d\d.txt',f)))]
    L = len(list_of_files)            
    
    data_spectrum = []
    name_spectrum = []
    specs = []
        
    print(L, 'calibration spectra were acquired.')
    
    for k in range(L):
        name = os.path.join(folder,list_of_files[k])
        data_spectrum = np.loadtxt(name)
        name_spectrum.append(list_of_files[k])
        specs.append(data_spectrum)
    
    wavelength_filepath = os.path.join(folder,wavelength_filename[0])
    londa = np.loadtxt(wavelength_filepath)
    
    start_notch = closer(londa, start_notch)
    end_notch = closer(londa, end_notch)
    start_power = closer(londa, start_power) 
    end_power = closer(londa, end_power)
    lower_londa_index = closer(londa, lower_londa)
    upper_londa_index = closer(londa, upper_londa)
    
    # ALLOCATING
    line_spec_raw = np.zeros((calibration_image_size_px,camera_px_length))
    line_spec_smooth = np.zeros((calibration_image_size_px,camera_px_length))
    line_spec = np.zeros((calibration_image_size_px,camera_px_length))
   
    for i in range(calibration_image_size_px):
        line_spec_raw[i,:] = np.array(specs[i])
    del specs
    
    ######################## SMOOTH ############################
    ######################## SMOOTH ############################
    ######################## SMOOTH ############################
    
    # SPLIT SIGNALS INTO STOKES AND ANTI-STOKES, longer range for smoothing
    line_stokes_raw = line_spec_raw[:,end_notch:]
    londa_stokes = londa[end_notch:]
    
    line_antistokes_raw = line_spec_raw[:,:start_notch]
    londa_antistokes = londa[:start_notch]
    
    # SMOOTHING
    print('Smoothing signals...')
    aux_line_stokes_smooth = sig.savgol_filter(line_stokes_raw, 
                                               window, deg, axis = 1, 
                                               mode=mode)
    aux_line_antistokes_smooth = sig.savgol_filter(line_antistokes_raw, 
                                               window, deg, axis = 1, 
                                               mode=mode)
    
    for i in range(repetitions-1):
        aux_line_stokes_smooth = sig.savgol_filter(aux_line_stokes_smooth,
                                                   window, deg, axis = 1, 
                                                   mode=mode)
        aux_line_antistokes_smooth = sig.savgol_filter(aux_line_antistokes_smooth,
                                                   window, deg, axis = 1, 
                                                   mode=mode)
    # Merge
    line_stokes_smooth = aux_line_stokes_smooth
    line_antistokes_smooth = aux_line_antistokes_smooth
    line_spec_smooth[:,end_notch:] = line_stokes_smooth
    line_spec_smooth[:,:start_notch] = line_antistokes_smooth
    
    print('Saving single plots of the calibration spectra as measured...')
    for i in range(calibration_image_size_px):
        plt.figure()
        pixel_name = 'i%02d_j00' % (i)
        plt.plot(londa, line_spec_raw[i,:], color='C0', linestyle='-', label='As measured')
        plt.plot(londa_stokes, line_stokes_smooth[i,:], color='k', linestyle='-', label='Smoothed')
        plt.plot(londa_antistokes, line_antistokes_smooth[i,:], color='k', linestyle='-')
        plt.legend()
        ax = plt.gca()
        ax.set_xlabel(r'Wavelength (nm)')
        ax.set_ylabel('Intensity (a.u.)')
        ax.axvline(londa[lower_londa_index], ymin = 0, ymax = 1, color='k', linestyle='--')
        ax.axvline(londa[upper_londa_index], ymin = 0, ymax = 1, color='k', linestyle='--')
        ax.axvline(londa[start_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
        ax.axvline(londa[end_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
        ax.axvline(londa[start_power], ymin = 0, ymax = 1, color='k', linestyle='--')
        ax.axvline(londa[end_power], ymin = 0, ymax = 1, color='k', linestyle='--')
        aux_folder = manage_save_directory(save_folder,'calibration_spectra_as_measured')
        figure_name = os.path.join(aux_folder, 'calibration_spec_%s_%s.png' % (pixel_name, NP))
        ax.set_ylim([y_min_spec_cal, y_max_spec_cal])
        plt.savefig(figure_name)
        plt.close()
    
    ######################## KILL NOTCH RANGE ############################
    ######################## KILL NOTCH RANGE ############################
    ######################## KILL NOTCH RANGE ############################

    line_spec_smooth[:,start_notch:end_notch] = np.nan
    
    ######################## FIND BKG AND MAX ############################
    ######################## FIND BKG AND MAX ############################
    ######################## FIND BKG AND MAX ############################
    
    # LOOK FOR MAX AND MIN IN STOKES RANGE
    line_stokes_smooth = line_spec_smooth[:,start_power:end_power]

    aux_sum = np.sum(line_stokes_smooth, axis=1)

    print('Finding max and bkg (min) spectra of the calibration...')            
    
    imin = np.argmin(aux_sum)
    bkg_smooth = line_spec_smooth[imin, :]
    bkg_raw = line_spec_raw[imin, :]
    noise_rms = np.sqrt(np.nansum(bkg_raw**2)/len(bkg_raw))
    
    imax = np.argmax(aux_sum)
    max_smooth = line_spec_smooth[imax, :]
    max_raw = line_spec_raw[imax, :]
    signal_rms = np.sqrt(np.nansum(max_smooth**2)/len(max_smooth))
    
    signal_to_background_ratio = (signal_rms/noise_rms)**2
    print('Signal RMS (max) to bkg RMS (min) ratio:', signal_to_background_ratio)
                
    # BACKGROUND SPECTRUM
    plt.figure()
    plt.plot(londa, bkg_raw)
    plt.plot(londa, bkg_smooth, '-k')
    ax = plt.gca()
    ax.set_xlabel(r'Wavelength (nm)')
    ax.set_ylabel('Intensity (a.u.)')
    ax.set_ylim([y_min_spec_cal, y_max_spec_cal])
    ax.axvline(londa[start_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[end_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[start_power], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[end_power], ymin = 0, ymax = 1, color='k', linestyle='--')
    figure_name = os.path.join(save_folder,'calibration_bkg_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'calibration/bkg')
    figure_name = os.path.join(aux_folder,'calibration_bkg_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()
    
    # MAX SPECTRUM
    plt.figure()
    plt.plot(londa, (max_raw - y_min_spec_cal)/10000, label = 'Raw')
    plt.plot(londa, (max_smooth - y_min_spec_cal)/10000, '-k', label = 'Smooth')
    ax = plt.gca()
    ax.legend(loc='upper left', prop={'size':13})
    ax.set_xlabel(r'Wavelength (nm)', fontsize=20)
    ax.set_ylabel('Photoluminescence (a.u.)', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=18)
#    ax.set_ylim([y_min_spec_cal, y_max_spec_cal])
    ax.set_ylim([0, (y_max_spec_cal - y_min_spec_cal)/10000])
    ax.text(582, 16, 'G', fontsize=16)
#    ax.axvline(londa[start_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
#    ax.axvline(londa[end_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
#    ax.axvline(londa[start_power], ymin = 0, ymax = 1, color='k', linestyle='--')
#    ax.axvline(londa[end_power], ymin = 0, ymax = 1, color='k', linestyle='--')
    figure_name = os.path.join(save_folder,'calibration_max_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'calibration/max')
    figure_name = os.path.join(aux_folder,'calibration_max_%s.png' % NP)
    plt.savefig(figure_name, dpi=300, bbox_plot = 'tight')
    figure_name = os.path.join(aux_folder,'calibration_max_%s.pdf' % NP)
    plt.savefig(figure_name, dpi=300, bbox_plot = 'tight', format='pdf')
    plt.close()
    
    # BACKGROUND SUBSTRACTION
    line_spec = line_spec_smooth - bkg_smooth
    
    ######################## DETERMINE COUNTS PER IRRADIANCE ############################
    ######################## DETERMINE COUNTS PER IRRADIANCE ############################
    ######################## DETERMINE COUNTS PER IRRADIANCE ############################

    print('Calculating counts per unit of time...')                
    
    aux_sum_stokes_smooth_corrected = np.sum(line_spec[:,start_power:end_power], axis=1)

    integrated_cts = np.max(aux_sum_stokes_smooth_corrected)
    integrated_cts_per_sec = integrated_cts/exp_time_fast
    
    print('Stokes-integrated signal (bkg corrected): %.0f +/- %.0f cts/s' % (integrated_cts_per_sec, 
          np.sqrt(integrated_cts_per_sec)))

    meas_pow_sample = factor*meas_pow_bfp
    err_meas_pow_sample = instrumental_err*meas_pow_sample
    irrad_calc = 2*meas_pow_sample/(np.pi*w0**2)
    err_irrad_calc = np.sqrt( (2/(np.pi*w0**2)*err_meas_pow_sample)**2 +
                             (4*meas_pow_sample*err_w0/(np.pi*w0**3))**2 )
    
    print('Irradiance max (calc): %.2f +/- %.2f mW/um2\n' % (irrad_calc, err_irrad_calc))

    #################### SAVE DATA ############################
    #################### SAVE DATA ############################
    #################### SAVE DATA ############################
    
    print('Saving calibration data...')
    
    aux_folder = manage_save_directory(save_folder,'calibration')
    to_save = [integrated_cts_per_sec, irrad_calc, err_irrad_calc]
    path_to_save = os.path.join(aux_folder,'integrated_counts_%s.dat' % NP)
    np.savetxt(path_to_save, to_save, fmt='%.6e')
    
    aux_folder = manage_save_directory(save_folder,'calibration')
    to_save = bkg_smooth
    path_to_save = os.path.join(aux_folder,'calibration_bkg_%s.dat' % NP)
    np.savetxt(path_to_save, to_save, fmt='%.6e')
    
    aux_folder = manage_save_directory(save_folder,'calibration')
    to_save = max_smooth
    path_to_save = os.path.join(aux_folder,'calibration_max_%s.dat' % NP)
    np.savetxt(path_to_save, to_save, fmt='%.6e')
    
    return