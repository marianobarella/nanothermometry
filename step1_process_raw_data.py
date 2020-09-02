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
from matplotlib import cm
from matplotlib_scalebar.scalebar import ScaleBar
import scipy.signal as sig
from functions_for_photoluminiscence import manage_save_directory, \
    classification, calc_r2, closer, lorentz2, fit_lorentz2

try:
    plt.style.use('for_confocal.mplstyle')
except:
    print('Pre-defined matplotlib style was not loaded.')

plt.ioff()
plt.close('all')

def process_power_during_confocal(folder, path_to):
    
    NP = folder.split('Spectrum_')[-1]
    
    common_path = os.path.join(path_to,'common_plots')

    save_folder = os.path.join(path_to, NP)
    
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
            
    list_of_files = os.listdir(folder)
    monitor_filename = [f for f in list_of_files if re.search('Trace_BS',f)]

    power_at_bs_filepath = os.path.join(folder, monitor_filename[0])
    data = np.loadtxt(power_at_bs_filepath)
    time = data[:,0]
    signal = data[:,1]
    
    avg_signal = np.mean(signal)
    std_signal = np.std(signal, ddof = 1)
    
    plt.figure()
    plt.plot(time, signal, '.-', color='C0')
    ax = plt.gca()
    ax.set_xlabel(r'Time (s)')
    ax.set_ylabel('Intensity (a.u.)')
#    ax.set_ylim([0.140, 0.170])
    ax.axhline(avg_signal, xmin = time[0], xmax = time[-1], color='C3', linestyle='-')
    ax.axhline(avg_signal + std_signal, xmin = time[0], xmax = time[-1], color='k', linestyle='--')
    ax.axhline(avg_signal - std_signal, xmin = time[0], xmax = time[-1], color='k', linestyle='--')
    aux_folder = manage_save_directory(save_folder,'power_monitor')
    figure_name = os.path.join(aux_folder, 'power_monitor_vs_time_%s.png' % (NP))
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'power_monitor')
    figure_name = os.path.join(aux_folder, 'power_monitor_vs_time_%s.png' % (NP))
    plt.savefig(figure_name)
    plt.close()
        
    print('Power monitor data was processed.')  
    print('monitor avg intensity = %.3f V' % avg_signal)  
    print('Ratio power std/avg = %.3f' % (std_signal/avg_signal))
    
    aux_folder = manage_save_directory(save_folder,'power_monitor')
    to_save = [avg_signal, std_signal]
    path_to_save = os.path.join(aux_folder, 'power_monitor_%s.dat' % NP)
    np.savetxt(path_to_save, to_save, fmt='%.3e')

    return

def process_confocal_to_bins(folder, path_to, totalbins, image_size_px, 
                             image_size_um, camera_px_length, window, deg, 
                             repetitions, mode, factor, meas_pow_bfp, start_notch, 
                             end_notch, start_power, end_power, start_spr, max_spr, 
                             lower_londa, upper_londa, exp_time_slow, w0_input, 
                             use_calibration_flag, last_bin_is_bkg_flag, plot_flag=False):
    
    pixel_size = image_size_um/image_size_px # in um^2
    print('Pixel size of spectral confocal image = %.3f um' % pixel_size)
    
    # plot limits
    y_min_spec_conf = 8000
    y_max_spec_conf = 100000
    
    NP = folder.split('Spectrum_')[-1]
    
    save_folder = os.path.join(path_to, NP)
    
    calibration_path = os.path.join(save_folder,'calibration')
    
    calibration_irrrad_filepath = os.path.join(calibration_path, 'integrated_counts_%s.dat' % NP)
    [cal_integrated_cts_per_sec, 
     cal_irrad_calc,
     err_cal_irrad_calc] = np.loadtxt(calibration_irrrad_filepath)

    calibration_max_filepath = os.path.join(calibration_path,'calibration_max_%s.dat' % NP)
    calibration_max_smooth = np.loadtxt(calibration_max_filepath)

    calibration_bkg_filepath = os.path.join(calibration_path,'calibration_bkg_%s.dat' % NP)
    calibration_bkg_smooth = np.loadtxt(calibration_bkg_filepath)

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
        
    print(L, 'spectra were acquired.')
    
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
    matrix_spec_raw = np.zeros((image_size_px,image_size_px,camera_px_length))
    matrix_spec_smooth = np.zeros((image_size_px,image_size_px,camera_px_length))
    matrix_spec = np.zeros((image_size_px,image_size_px,camera_px_length))
    matrix_spec_normed = np.zeros((image_size_px,image_size_px,camera_px_length))
    binned_power = np.zeros((image_size_px,image_size_px))
   
    for i in range(image_size_px):
        for j in range(image_size_px):
            matrix_spec_raw[i,j,:] = np.array(specs[i*image_size_px+j])
    del specs
    
    ######################## SMOOTH ############################
    ######################## SMOOTH ############################
    ######################## SMOOTH ############################
    
    # SPLIT SIGNALS INTO STOKES AND ANTI-STOKES
    matrix_stokes_raw = matrix_spec_raw[:,:,end_notch:]
    londa_stokes = londa[end_notch:]
    
    matrix_antistokes_raw = matrix_spec_raw[:,:,:start_notch]
    londa_antistokes = londa[:start_notch]
    
    # SMOOTHING
    print('Smoothing signals...')
    aux_matrix_stokes_smooth = sig.savgol_filter(matrix_stokes_raw, 
                                               window, deg, axis = 2, 
                                               mode=mode)
    aux_matrix_antistokes_smooth = sig.savgol_filter(matrix_antistokes_raw, 
                                               window, deg, axis = 2, 
                                               mode=mode)
    
    for i in range(repetitions - 1):
        aux_matrix_stokes_smooth = sig.savgol_filter(aux_matrix_stokes_smooth,
                                                   window, deg, axis = 2, 
                                                   mode=mode)
        aux_matrix_antistokes_smooth = sig.savgol_filter(aux_matrix_antistokes_smooth,
                                                   window, deg, axis = 2, 
                                                   mode=mode)
    # Merge
    matrix_stokes_smooth = aux_matrix_stokes_smooth
    matrix_antistokes_smooth = aux_matrix_antistokes_smooth
    matrix_spec_smooth[:,:,end_notch:] = matrix_stokes_smooth
    matrix_spec_smooth[:,:,:start_notch] = matrix_antistokes_smooth
    
    if plot_flag:
        print('Saving single spectrums as measured...')
        for i in range(image_size_px):
            for j in range(image_size_px):
                plt.figure()
                pixel_name = 'i%02d_j%02d' % (i,j)
                plt.plot(londa, matrix_spec_raw[i,j,:], color='C0', linestyle='-', label='As measured')
                plt.plot(londa_stokes, matrix_stokes_smooth[i,j,:], color='k', linestyle='-', label='Smoothed')
                plt.plot(londa_antistokes, matrix_antistokes_smooth[i,j,:], color='k', linestyle='-')
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
                aux_folder = manage_save_directory(save_folder,'pl_spectra_as_measured')
                figure_name = os.path.join(aux_folder, 'spec_%s_%s.png' % (pixel_name, NP))
                ax.set_ylim([y_min_spec_conf,y_max_spec_conf])
                plt.savefig(figure_name)
                plt.close()
        print('Saving all spectra as measured in one figure.')
        aux_matrix_spec_smooth2 = sig.savgol_filter(matrix_spec_raw, 
                                               window, deg, axis = 2, 
                                               mode=mode)
        for i in range(repetitions - 1):
            aux_matrix_spec_smooth2 = sig.savgol_filter(aux_matrix_spec_smooth2,
                                                   window, deg, axis = 2, 
                                                   mode=mode)
        plt.figure()
        for i in range(image_size_px):
            for j in range(image_size_px):
                plt.plot(londa, aux_matrix_spec_smooth2[i,j,:])
#                plt.legend()
        ax = plt.gca()
#                ax.set_ylim([y_min_spec_conf,y_max_spec_conf])
        ax.set_xlabel(r'Wavelength (nm)', fontsize=20)
#                ax.set_ylabel('Intensity (bkg corrected) (a.u.)' , fontsize=20)
        ax.set_ylabel('Intensity (a.u.)' , fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=18)
        figure_name = os.path.join(save_folder,'all_spec_as_measured_%s.png' % (NP))
        plt.savefig(figure_name, dpi=300)
        figure_name = os.path.join(save_folder,'all_spec_as_measured_%s.pdf' % (NP))
        plt.savefig(figure_name, dpi=300, format='pdf')
        plt.close()  
    
    ######################## KILL NOTCH RANGE ############################
    ######################## KILL NOTCH RANGE ############################
    ######################## KILL NOTCH RANGE ############################

    matrix_spec_smooth[:,:,start_notch:end_notch] = np.nan
    
    ######################## FIND BKG AND MAX ############################
    ######################## FIND BKG AND MAX ############################
    ######################## FIND BKG AND MAX ############################
    
    # LOOK FOR MAX AND MIN IN STOKES RANGE
    matrix_stokes_smooth = matrix_spec_smooth[:,:,start_power:end_power]
    aux_sum = np.sum(matrix_stokes_smooth, axis=2)
    
    print('Finding max and bkg (min) spectra...')            
    
    imin, jmin = np.unravel_index(np.argmin(aux_sum, axis=None), aux_sum.shape)
    bkg_smooth = matrix_spec_smooth[imin, jmin, :]
    bkg_raw = matrix_spec_raw[imin, jmin, :]
    noise_rms = np.sqrt(np.nansum(bkg_raw**2)/len(bkg_raw))
    
    imax, jmax = np.unravel_index(np.argmax(aux_sum, axis=None), aux_sum.shape)
    max_smooth = matrix_spec_smooth[imax, jmax, :]
    max_raw = matrix_spec_raw[imax, jmax, :]
    signal_rms = np.sqrt(np.nansum(max_smooth**2)/len(max_smooth))
    
    signal_to_background_ratio = (signal_rms/noise_rms)**2
    print('Signal RMS (max) to bkg RMS (min) ratio:', signal_to_background_ratio)
    
    integrated_calibration_bkg_stokes_smooth = np.sum(calibration_bkg_smooth[start_power:end_power])
    integrated_bkg_stokes_smooth = np.sum(bkg_smooth[start_power:end_power])
    bkg_to_bkg_ratio = integrated_bkg_stokes_smooth/integrated_calibration_bkg_stokes_smooth
    print('Bkg-confocal to bkg-calibration ratio:', bkg_to_bkg_ratio)
    print('Warning! If fast and slow exposure times differ, this ratio would be far from 1. \n Base line clamp EMCDD function would cause this.')
    
    signal_to_background = max_raw/bkg_raw
    # signal_to_background SPECTRUM
    plt.figure()
    plt.plot(londa, signal_to_background)
    ax = plt.gca()
    ax.set_xlabel(r'Wavelength (nm)')
    ax.set_ylabel('Intensity (a.u.)')
    ax.axvline(londa[lower_londa_index], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[upper_londa_index], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[start_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[end_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[start_power], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[end_power], ymin = 0, ymax = 1, color='k', linestyle='--')
    figure_name = os.path.join(save_folder,'signal_to_background_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()
                
    # BACKGROUND SPECTRUM
    plt.figure()
    plt.plot(londa, bkg_raw, label='Raw')
    plt.plot(londa, bkg_smooth, '-k', label='Smooth')
    plt.plot(londa, calibration_bkg_smooth, ':', color = 'C3', label='Bkg calib.')
    ax = plt.gca()
    ax.set_xlabel(r'Wavelength (nm)')
    ax.set_ylabel('Intensity (a.u.)')
    ax.set_ylim([y_min_spec_conf,y_max_spec_conf])
    ax.axvline(londa[lower_londa_index], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[upper_londa_index], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[start_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[end_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[start_power], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[end_power], ymin = 0, ymax = 1, color='k', linestyle='--')
    plt.legend(loc='upper left', prop = {'size':8})
    figure_name = os.path.join(save_folder,'bkg_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'bkg')
    figure_name = os.path.join(aux_folder,'bkg_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()
    
    # MAX SPECTRUM
    plt.figure()
    plt.plot(londa, max_raw, label='Raw')
    plt.plot(londa, max_smooth, '-k', label='Smooth')
    plt.plot(londa, calibration_max_smooth, ':', color = 'C3', label='Max calib.')
    ax = plt.gca()
    ax.set_xlabel(r'Wavelength (nm)')
    ax.set_ylabel('Intensity (a.u.)')
    ax.set_ylim([y_min_spec_conf,y_max_spec_conf])
    ax.axvline(londa[lower_londa_index], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[upper_londa_index], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[start_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[end_notch], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[start_power], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[end_power], ymin = 0, ymax = 1, color='k', linestyle='--')
    plt.legend(loc='upper left', prop = {'size':8})
    figure_name = os.path.join(save_folder,'max_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'max')
    figure_name = os.path.join(aux_folder,'max_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()

    # BACKGROUND SUBSTRACTION AND NORMALIZARION
    # BACKGROUND SUBSTRACTION AND NORMALIZARION
    # BACKGROUND SUBSTRACTION AND NORMALIZARION

    if use_calibration_flag:
        print('Background substraction and normalization (for Stokes check) using calibration...')
        matrix_spec = matrix_spec_smooth - calibration_bkg_smooth
        matrix_spec_normed = (matrix_spec_smooth - calibration_bkg_smooth) / (calibration_max_smooth - calibration_bkg_smooth)
    else:
        print('Background substraction and normalization (for Stokes check) without using calibration. Using the confocal itself...')
        matrix_spec = matrix_spec_smooth - bkg_smooth
        matrix_spec_normed = (matrix_spec_smooth - bkg_smooth) / (max_smooth - bkg_smooth)

    if plot_flag:
        print('Saving single spectrums without bkg...')
        for i in range(image_size_px):
            for j in range(image_size_px):
                plt.figure()
                pixel_name = 'i%02d_j%02d' % (i,j)
                plt.plot(londa, matrix_spec[i,j,:], label=pixel_name)
#                plt.legend()
                ax = plt.gca()
#                ax.set_ylim([y_min_spec_conf,y_max_spec_conf])
                ax.set_xlabel(r'Wavelength (nm)', fontsize=20)
#                ax.set_ylabel('Intensity (bkg corrected) (a.u.)' , fontsize=20)
                ax.set_ylabel('Intensity (a.u.)' , fontsize=20)
                plt.tick_params(axis='both', which='major', labelsize=18)
                aux_folder = manage_save_directory(save_folder,'pl_spectra_minus_bkg')
                figure_name = os.path.join(aux_folder,'spec_minus_bkg_%s_%s.png' % (pixel_name, NP))
                plt.savefig(figure_name, dpi=300)
#                figure_name = os.path.join(aux_folder,'spec_minus_bkg_%s_%s.pdf' % (pixel_name, NP))
#                plt.savefig(figure_name, dpi=300, format='pdf')
                plt.close()  

    print('Plotting stokes spectral quotient to CHECK for proportionality...')
    plt.figure()
    for i in range(image_size_px):
        for j in range(image_size_px):
            pixel_name = 'i%02d_j%02d' % (i,j)
            plt.plot(londa, matrix_spec_normed[i,j,:], label=pixel_name)
    ax = plt.gca()
    ax.set_xlabel(r'Wavelength (nm)')
    ax.set_ylabel('Normalized Intensity (a.u.)')
    ax.axvline(londa[end_notch], ymin = -1, ymax = 2, color ='k', linestyle='--')
    ax.axvline(londa[start_power], ymin = 0, ymax = 1, color='k', linestyle='--')
    ax.axvline(londa[end_power], ymin = -1, ymax = 2, color ='k', linestyle='--')
    ax.set_ylim([-0.1,1.1])
    ax.set_xlim([0.99*londa[end_notch], 1.01*londa[-1]])
    aux_folder = manage_save_directory(save_folder,'pl_stokes_normalized')
    figure_name = os.path.join(aux_folder,'spec_stokes_normalized_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'pl_stokes_normalized')
    figure_name = os.path.join(aux_folder,'spec_stokes_normalized_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()

    ######################## FITTING LORENTZ CURVE AT STOKES ############################
    ######################## FITTING LORENTZ CURVE AT STOKES ############################
    ######################## FITTING LORENTZ CURVE AT STOKES ############################

    print('Extracting SPR peak and width...')
    londa_spr = londa[start_spr:]
    sum_stokes_spr = np.sum(matrix_spec[:,:,start_spr:], axis = (0, 1))
    sum_stokes_spr = sum_stokes_spr/np.max(sum_stokes_spr)
    # Fitting the data
    # initial parameter guesses
    # [amplitude, FWHM, center, offset]
    try:        
        init_params = np.array([1, 50, max_spr, 0.05], dtype=np.double)
        # Get the fitting parameters for the best lorentzian
        best_lorentz, err = fit_lorentz2(init_params, londa_spr, sum_stokes_spr)
        # calculate the errors
        lorentz_fitted = lorentz2(londa_spr, *best_lorentz)
        r2_coef_pearson = calc_r2(sum_stokes_spr, lorentz_fitted)
        full_lorentz_fitted = lorentz2(londa, *best_lorentz)
        londa_max_pl = best_lorentz[2]
        width_pl = best_lorentz[1]
        print('SPR wavelenth (max) = %.2f nm, Width = %.2f nm' % (londa_max_pl, width_pl))
    except RuntimeError:
        print('SPR fitting did not converge. Analysis must go on... ignoring NP\'s SPR.')
        r2_coef_pearson = 999
        full_lorentz_fitted = np.zeros(len(londa))
        londa_max_pl = 0
        width_pl = 0
    
    # STOKES FITTING
    plt.figure()
    plt.plot(londa_spr, sum_stokes_spr,'C0', label='Data %s' % NP)
    plt.plot(londa, full_lorentz_fitted,'k--', label='Lorentz fit')
    plt.legend()
    ax = plt.gca()
    ax.set_xlabel(r'Wavelength (nm)')
    ax.set_ylabel('Sum of intensities (a.u.)')
    plt.ylim([0,1.05])
    aux_folder = manage_save_directory(save_folder,'spr')
    figure_name = os.path.join(aux_folder, 'sum_specs_mean_fitted_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'spr')
    figure_name = os.path.join(aux_folder, 'sum_specs_mean_fitted_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close() 

    ######################## IRRADIANCE ############################
    ######################## IRRADIANCE ############################
    ######################## IRRADIANCE ############################

    print('Irradiance (calc): %.2f +/- %.2f mW/um2' % \
          (cal_irrad_calc, err_cal_irrad_calc))
    print('Calibration returned (max, bkg corrected): %.0f cts/s' % cal_integrated_cts_per_sec)

    matrix_sum_stokes_spec = np.sum(matrix_spec[:,:,start_power:end_power], axis=2)
    matrix_integrated_cts_per_sec = matrix_sum_stokes_spec/exp_time_slow
    
    imax, jmax = np.unravel_index(np.argmax(matrix_integrated_cts_per_sec, axis=None), matrix_integrated_cts_per_sec.shape)
    imin, jmin = np.unravel_index(np.argmin(matrix_integrated_cts_per_sec, axis=None), matrix_integrated_cts_per_sec.shape)
    
    max_integrated_cts_per_sec = matrix_integrated_cts_per_sec[imax, jmax]
    min_integrated_cts_per_sec = matrix_integrated_cts_per_sec[imin, jmin]
    
    # Poisson
    err_max_integrated_cts_per_sec = np.sqrt(max_integrated_cts_per_sec)
    err_min_integrated_cts_per_sec = np.sqrt(min_integrated_cts_per_sec)
    
    print('Confocal image returned (max, bkg corrected): %.0f +/- %.0f cts/s' % \
          (max_integrated_cts_per_sec, err_max_integrated_cts_per_sec))
    print('Confocal image returned (min, bkg corrected): %.0f +/- %.0f cts/s' % \
          (min_integrated_cts_per_sec, err_min_integrated_cts_per_sec))
    
    ratio_confocal_max_to_calibration_max = max_integrated_cts_per_sec/cal_integrated_cts_per_sec
    print('Ratio max slow/max fast:', ratio_confocal_max_to_calibration_max)
    
    if use_calibration_flag:
        print('Calculating irradiance per pixel using calibration...')
        pixel_cts_per_sec = matrix_integrated_cts_per_sec/cal_integrated_cts_per_sec
    else:
        print('Calculating irradiance per pixel without calibration. Using the confocal itself...')
        pixel_cts_per_sec = matrix_integrated_cts_per_sec/max_integrated_cts_per_sec
    
    
    pixel_irrad_estimated = pixel_cts_per_sec*cal_irrad_calc
    err_pixel_irrad_estimated = pixel_cts_per_sec * np.abs(err_cal_irrad_calc)
    
    # IMAGE POWER HOT MAP USING CALIBRATION
    plt.figure()
    img = plt.imshow(pixel_irrad_estimated, interpolation = 'none', cmap = 'magma')#'inferno')#'viridis')#'hot')
    ax = plt.gca()
    ax.set_xticks(range(image_size_px))
    ax.set_yticks(range(image_size_px))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    scalebar = ScaleBar(pixel_size*1e-6, location = 'upper left') 
    ax.add_artist(scalebar)
    cbar = plt.colorbar()
    cbar.ax.set_title(u'Irradiance (mW/µm$^{2}$)', fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    figure_name = os.path.join(save_folder,'pixel_irrad_estimated_from_cal_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300)
    aux_folder = manage_save_directory(common_path,'pixel_irrad')
    figure_name = os.path.join(aux_folder,'pixel_irrad_estimated_from_cal_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300)
    figure_name = os.path.join(aux_folder,'pixel_irrad_estimated_from_cal_%s.pdf' % NP)
    plt.savefig(figure_name, dpi = 300, format='pdf')
    plt.close()

    ######################## BINNING ############################
    ######################## BINNING ############################
    ######################## BINNING ############################
    
    print('Binning...')

    for i in range(image_size_px):
        for j in range(image_size_px):
            nbin = classification(pixel_cts_per_sec[i,j], totalbins, [0,1])
            binned_power[i,j] = nbin
    
    plt.figure()
    ax = plt.gca()
    ax.set_xticks(range(image_size_px))
    ax.set_yticks(range(image_size_px))
    ax.set_xticks([])
    ax.set_yticks([])
    rango = range(0, totalbins + 1)
    ticks_label_rango = [float(i)/totalbins for i in rango]
    img = ax.imshow(binned_power, interpolation = 'none', cmap = cm.plasma)
    cbar = plt.colorbar(img, ax = ax, boundaries = rango, ticks = rango)
    cbar.ax.set_yticklabels(ticks_label_rango, fontsize = 20)  # vertically oriented colorbar
    plt.grid(False)
    figure_name = os.path.join(save_folder,'binned_power_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300)
    aux_folder = manage_save_directory(common_path,'binned_power')
    figure_name = os.path.join(aux_folder,'binned_power_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300)
    figure_name = os.path.join(aux_folder,'binned_power_%s.pdf' % NP)
    plt.savefig(figure_name, dpi = 300, format = 'pdf')
    plt.close()
    
    # SPECTRUMS IN SAME GRAPH ACCORDING TO ITS POWER BIN
    mean_specs = np.zeros((totalbins, camera_px_length))
    mean_irrad = np.zeros(totalbins)
    err_mean_irrad = np.zeros(totalbins)
    hist_bin = np.zeros(totalbins)
    for s in range(totalbins):
        u, v = np.where(binned_power == s)
        counter = 0
        aux_spec = np.zeros(camera_px_length)
        aux_irrad = 0
        aux_err_irrad = 0
        plt.figure()
        for i, j in zip(u, v):
            plt.plot(londa, matrix_spec[i,j,:])
            aux_irrad += pixel_irrad_estimated[i,j]
            aux_err_irrad += err_pixel_irrad_estimated[i,j]**2
            aux_spec += matrix_spec[i,j,:]
            counter += 1
        if counter == 0:
            print('Bin %d has NO spectra to average.' % s)
            aux_spec = aux_spec
            aux_irrad = aux_irrad
            aux_err_irrad = aux_err_irrad
        else:
            print('Bin %d has %d spectra to average.' % (s, counter))
            aux_spec = aux_spec/counter
            aux_irrad = aux_irrad/counter
            aux_err_irrad = np.sqrt(aux_err_irrad)/counter
        mean_specs[s] = aux_spec
        mean_irrad[s] = aux_irrad
        err_mean_irrad[s] = aux_err_irrad
        hist_bin[s] = counter
        plt.plot(londa, aux_spec, '--k')
        plt.title('Bin %d' % s)
        plt.ylim([0, y_max_spec_conf - y_min_spec_conf])
        aux_folder = manage_save_directory(save_folder, 'pl_in_bins')
        figure_name = os.path.join(aux_folder, 'bin_%d_%s.png' % (s, NP))
        plt.savefig(figure_name)
        plt.close()
    
    # binning error determination
    estimated_value = hist_bin/np.sum(hist_bin)
    err_bin = np.round(np.sqrt(hist_bin*(1 - estimated_value))) # binomial std deviation
    
    # find theorical binning
    edges = np.array(range(totalbins))
    centers = edges + 0.5
    centers_distro = centers[::-1]
    height_distro = np.zeros(totalbins)
    step = 1/totalbins
    binned_range = np.arange(0, 1 + step, step)[::-1]
    for i in range(totalbins-1):
        I_1 = binned_range[i]
        I_2 = binned_range[i+1]
        height_distro[i] = 0.5*(w0_input**2)*np.log(I_1/I_2)*np.pi/(pixel_size**2)
#        print(height_distro[i], 'between', I_1, I_2)
    height_distro[-1] = image_size_px**2 - np.sum(height_distro)
    
    # plot histogram of binning
    plt.figure()
    plt.bar(centers, hist_bin, align='center')
    plt.errorbar(centers, hist_bin, fmt='', linestyle='', yerr=err_bin)
    plt.plot(centers_distro, height_distro, 'o-', color='C3')
    plt.xlabel('Bin')
    plt.ylabel('Number of spectra')
    plt.ylim([0,image_size_px**2])
    aux_folder = manage_save_directory(save_folder,'pl_in_bins')
    figure_name = os.path.join(aux_folder,'hist_of_binning_%s_lin.png' % (NP))
    plt.savefig(figure_name)
    plt.ylim([0.9,image_size_px**2])
    ax = plt.gca()
    ax.set_yscale('log')
    aux_folder = manage_save_directory(save_folder,'pl_in_bins')
    figure_name = os.path.join(aux_folder,'hist_of_binning_%s_log.png' % (NP))
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'hist_of_binning')
    figure_name = os.path.join(aux_folder,'hist_of_binning_%s_log.png' % (NP))
    plt.savefig(figure_name)
    plt.close()

    ################### LAST BIN IS BKG ########################
    ################### LAST BIN IS BKG ########################
    ################### LAST BIN IS BKG ########################
    
    if last_bin_is_bkg_flag:        
        # correct signals to account for mean bkg (last bin)
        print('Correction to bin spectra is being applied...')
        print('Considering last bin as a bkg correction...')
        print('Bin 0 irradiation = %.3f' % mean_irrad[0])
        corrected_mean_specs = mean_specs - mean_specs[0]
        corrected_mean_irrad = mean_irrad - mean_irrad[0]
    else:
        corrected_mean_specs = mean_specs
        corrected_mean_irrad = mean_irrad

    max_max = np.nanmax(corrected_mean_specs, axis = (0,1))
#    print(max_max)
    
    rango = list(range(0, totalbins + 1))
    norm = plt.Normalize()
    colormap = plt.cm.plasma(norm(rango))

    plt.figure()
    for i in range(totalbins):
        color_iter = colormap[i]
        plt.plot(londa, corrected_mean_specs[i]/max_max, 
                 color = color_iter,
                 label='%.2f' % mean_irrad[i])
#                 label='Group %d' % i)
#       plt.plot(londa, corrected_mean_specs[i], label='Bin %d' % i)
    ax = plt.gca()
    ax.set_xlabel(r'Wavelength (nm)', fontsize=20)
    ax.set_ylabel('Intensity (a.u.)', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=18)
    plt.ylim([0, 1.02])
#    plt.ylim([0, y_max_spec_conf - y_min_spec_conf])
    plt.legend(loc='upper right', prop={'size':11})
    aux_folder = manage_save_directory(save_folder,'pl_in_bins')
    figure_name = os.path.join(aux_folder,'all_bins_corrected_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300)
    aux_folder = manage_save_directory(common_path,'all_bins')
    figure_name = os.path.join(aux_folder,'all_bins_corrected_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300)
    figure_name = os.path.join(aux_folder,'all_bins_corrected_%s.pdf' % NP)
    plt.savefig(figure_name, dpi = 300, format = 'pdf')
    plt.close()
    
    plt.figure()
    for i in range(totalbins):
        color_iter = colormap[i]
        plt.plot(londa, corrected_mean_specs[i]/max_max, 
                 color = color_iter,
                 label='Group %d' % i)
    ax = plt.gca()
    ax.set_xlabel(r'Wavelength (nm)', fontsize=20)
    ax.set_ylabel('Intensity (a.u.)', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=18)
#    plt.ylim([0, 0.2*(y_max_spec_conf - y_min_spec_conf)])
    plt.ylim([0, 0.225])
    plt.xlim([497, 523])
    plt.legend(loc='upper left', prop={'size':10})
#    ax.axvline(londa[lower_londa_index], ymin=0, ymax=1, linestyle='--', color='k')
#    ax.axvline(londa[upper_londa_index], ymin=0, ymax=1, linestyle='--', color='k')
    aux_folder = manage_save_directory(save_folder, 'antistokes_in_bins')
    figure_name = os.path.join(aux_folder,'all_bins_corrected_antistokes_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300)
    aux_folder = manage_save_directory(common_path, 'antistokes_in_bins')
    figure_name = os.path.join(aux_folder,'all_bins_corrected_antistokes_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300)
    figure_name = os.path.join(aux_folder,'all_bins_corrected_antistokes_%s.pdf' % NP)
    plt.savefig(figure_name, dpi = 300, format = 'pdf')
    plt.close()
    
    # FIGURE FOR PAPER, NORMALIZATION USING MAX STOKES PER BIN
    # FIGURE FOR PAPER, NORMALIZATION USING MAX STOKES PER BIN
    # FIGURE FOR PAPER, NORMALIZATION USING MAX STOKES PER BIN
    
    plt.figure()
    for i in range(totalbins):
        max_bin = np.nanmax(corrected_mean_specs[i])
        color_iter = colormap[i]
        plt.plot(londa, corrected_mean_specs[i]/max_bin, 
                 color = color_iter,
                 label='%.2f' % mean_irrad[i])
#                 label='Group %d' % i)
    ax = plt.gca()
    ax.set_xlabel(r'Wavelength (nm)', fontsize=20)
    ax.set_ylabel('Normalized Anti-Stokes (a.u.)', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=18)
    plt.xlim([505, 523])
    plt.legend(loc='lower right', prop={'size':12})
    ax.set_yscale('log')
    plt.ylim([0.005, 0.3])
    aux_folder = manage_save_directory(common_path, 'antistokes_in_bins')
    figure_name = os.path.join(aux_folder,'all_bins_corrected_nomalized_antistokes_log_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300)
    figure_name = os.path.join(aux_folder,'all_bins_corrected_nomalized_antistokes_log_%s.pdf' % NP)
    plt.savefig(figure_name, dpi = 300, format = 'pdf')
    plt.close()
    
    plt.figure()
    for i in range(totalbins):
        max_bin = np.nanmax(corrected_mean_specs[i])
        color_iter = colormap[i]
        plt.plot(londa, corrected_mean_specs[i]/max_bin, 
                 color = color_iter,
                 label='%.2f' % mean_irrad[i])
#                 label='Group %d' % i)
    ax = plt.gca()
    ax.set_xlabel(r'Wavelength (nm)', fontsize=20)
    ax.set_ylabel('Normalized intensity (a.u.)', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=18)
#    plt.xlim([497, 523])
    plt.legend(loc='lower center', prop={'size':11})
    ax.set_yscale('log')
    plt.ylim([0.004, 1.25])
    aux_folder = manage_save_directory(common_path, 'antistokes_in_bins')
    figure_name = os.path.join(aux_folder,'all_bins_corrected_nomalized_intensity_log1_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300)
    figure_name = os.path.join(aux_folder,'all_bins_corrected_nomalized_intensity_log1_%s.pdf' % NP)
    plt.savefig(figure_name, dpi = 300, format = 'pdf')
    plt.ylim([0.1, 1.05])
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('Normalized intensity (a.u.)', fontsize=20, rotation=270, labelpad=25)
    aux_folder = manage_save_directory(common_path, 'antistokes_in_bins')
    figure_name = os.path.join(aux_folder,'all_bins_corrected_nomalized_intensity_log2_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300)
    figure_name = os.path.join(aux_folder,'all_bins_corrected_nomalized_intensity_log2_%s.pdf' % NP)
    plt.savefig(figure_name, dpi = 300, format = 'pdf')
    plt.close()
    
    
    #################### SAVE DATA ############################
    #################### SAVE DATA ############################
    #################### SAVE DATA ############################
    print('Saving processed data...')
    
    aux_folder = manage_save_directory(save_folder,'pl_in_bins')
    to_save = corrected_mean_specs.T
    path_to_save = os.path.join(aux_folder,'all_bins_%s.dat' % NP)
    np.savetxt(path_to_save, to_save, fmt='%.3e')
    
    to_save = np.vstack((edges,
                         height_distro[::-1],
                         hist_bin)).T
    path_to_save = os.path.join(aux_folder,'bin_distro_%s.dat' % NP)
    np.savetxt(path_to_save, to_save, fmt='%.0f')
        
    to_save = londa
    path_to_save = os.path.join(aux_folder,'londa_%s.dat' % NP)
    np.savetxt(path_to_save, to_save, fmt='%.3e')
        
    to_save = np.array([corrected_mean_irrad, err_mean_irrad]).T
    path_to_save = os.path.join(aux_folder,'bin_irradiance_%s.dat' % NP)
    np.savetxt(path_to_save, to_save, fmt='%.3e')
    
    aux_folder = manage_save_directory(save_folder,'spr')
    to_save = [londa_max_pl, width_pl, r2_coef_pearson]
    path_to_save = os.path.join(aux_folder,'spr_fitted_parameters_%s.dat' % NP)
    np.savetxt(path_to_save, to_save, fmt='%.3e')
    
    return 