# -*- coding: utf-8 -*-
"""
Analysis of single AuNPs photoluminiscence spectra for temperature calculation 
acquired with PySpectrum at CIBION

Mariano Barella

14 jan 2020

"""

import os
import re
from functions_for_photoluminiscence import save_parameters
import step0_calibrate_irrad as step0
import step1_process_raw_data as step1
import step2_fit_antistokes_ratio as step2
import step3_calculate_temp as step3
import step4_statistics as step4

def run_analysis(T_lab_guess, beta_guess, \
                 find_beta, find_Tzero, \
                 meas_pow_bfp, NP_folder):
    # Parameters to load
    totalbins = 10 #number of bins
    zero_in_kelvin = 273 # in K
    Tzero_guess = zero_in_kelvin + T_lab_guess
    window, deg, repetitions = 51, 1, 2
    mode = 'interp'
    factor = 0.47 # factor de potencia en la muestra
    image_size_px = 12 # IN PIXELS
    image_size_um = 0.8 # IN um
    camera_px_length = 1002 # number of pixels (number of points per spectrum)
    calibration_image_size_px = 2 # IN PIXELS
    calibration_image_size_um = 0.8 # IN um
#    exp_time_fast = 10 # in seconds
    exp_time_fast = 10 # in seconds
    exp_time_slow = 2 # in seconds
    w0 = 0.342 # in um
    radius = 40 # in nm
    sigma_abs = 18498 # in nm^2 
    #radius = 76 # in nm
    #sigma_abs = 34138 # in nm^2 
#    radius = 51.5 # in nm
#    sigma_abs = 20949 # in nm^2 
#    radius = 27 # in nm
#    sigma_abs = 8599 # in nm^2 
    
    start_notch = 522 # in nm where notch starts ~525 nm (safe zone)
    end_notch = 543 # in nmwhere notch starts ~540 nm (safe zone)
    start_power = 550 # in nm from here we are going to calculate irradiance
    end_power = 570 # in nm from start_power up to this value we calculate irradiance
    start_spr = end_notch # in nm lambda from where to fit lorentz
    max_spr = 560 # wavelength of the spr maximum, initial parameter for the fitting algorithm 
    
    plot_flag = False # if True will save all spectra's plots for each pixel

    # Limites for the antistokes range
    lower_londa = 510 # 503 nm
    upper_londa = start_notch - 2
        
    # Parameters to load
    # Threhsold: check if data-point is going to be included i  n the analysis
    # If both crit0eria do not apply, erase datapoint
    alpha = 0.05 # alpha level for chi-squared test (compare with p-value), coarse criteria
    R2th = 0.8 # correlation coefficient threhsold, fine criteria
      
    single_NP_flag = False # if True only one NP will be analyzed
    NP_int = 3 # NP to be analyzed
    
    monitor_flag = False # flag to account for power monitor, False = not measured
    
    use_calibration_flag = False # flag to determine if the calibration (cts/s) are going to be used
    
    last_bin_is_bkg_flag = True # True if last bin is considered a background correction
        
    do_step0, do_step1, do_step2, do_step3, do_step4 = 1,1,1,1,1
#    do_step0, do_step1, do_step2, do_step3, do_step4 = 0,0,1,1,1
#    do_step0, do_step1, do_step2, do_step3, do_step4 = 0,1,0,0,0
#    do_step0, do_step1, do_step2, do_step3, do_step4 = 0,0,1,0,0
#     do_step0, do_step1, do_step2, do_step3, do_step4 = 0,0,0,0,1
#    do_step0, do_step1, do_step2, do_step3, do_step4 = 1,0,0,0,0
    
    base_folder = '/home/mariano/datos_mariano/posdoc/experimentos_PL_arg'
#    base_folder = '/home/marian/datos_vostro/posdoc/experimentos_PL_arg'
    
    parent_folder = os.path.join(base_folder, NP_folder)
    
    list_of_folders = os.listdir(parent_folder)
    list_of_folders = [f for f in list_of_folders if os.path.isdir(os.path.join(parent_folder,f))]
    list_of_folders_slow = [f for f in list_of_folders if re.search('Slow_Confocal_Spectrum',f)]
    list_of_folders_slow.sort()
    
    list_of_folders_fast = [f for f in list_of_folders if re.search('Fast_Confocal_Spectrum',f)]
    list_of_folders_fast.sort()
    
    path_to = os.path.join(parent_folder,'processed_data')
                
    save_parameters(path_to, totalbins, meas_pow_bfp, Tzero_guess, beta_guess, window, deg, repetitions ,\
                    mode, factor, image_size_px, image_size_um, camera_px_length, \
                    start_notch, end_notch, start_power, end_power, start_spr, max_spr, plot_flag, \
                    lower_londa, upper_londa, alpha, R2th, single_NP_flag, \
                    NP_int, do_step0, do_step1, do_step2, do_step3, do_step4)
    
    if find_beta and not find_Tzero:
        print('\nTzero is known. Finding beta...')
    elif not find_beta and find_Tzero:
        print('\nBeta is known. Finding Tzero...')
    
    if single_NP_flag:
        print('\nAnalyzing only NP %d.' % NP_int)
    
    #####################################################################
    #####################################################################
    #####################################################################
        
    if do_step0:
        for f in list_of_folders_fast:
            if single_NP_flag:    
                if not re.search('_NP_%03d' % NP_int, f):
                    continue
                else:
                    print('\nTarget NP found!')
                
            folder = os.path.join(parent_folder,f)
        
            print('\n === NP folder:',f)
    
            print('\n------------------------  STEP 0')
        
            step0.process_irrad_calibration(folder, path_to, calibration_image_size_px, 
                                           calibration_image_size_um, camera_px_length, window, 
                                           deg, repetitions, mode, factor, meas_pow_bfp, 
                                           start_notch, end_notch, start_power, end_power, 
                                           start_spr, lower_londa, upper_londa, exp_time_fast, w0)
    else:
        print('\nSTEP 0 was not executed.')
    
    if do_step1:   
        for f in list_of_folders_slow:
            if single_NP_flag:    
                if not re.search('_NP_%03d' % NP_int, f):
                    continue
                else:
                    print('\nTarget NP found!')
                
            folder = os.path.join(parent_folder,f)
        
            print('\n === NP folder:',f)
    
            print('\n------------------------  STEP 1')
            
            if monitor_flag:
                step1.process_power_during_confocal(folder, path_to)
        
            step1.process_confocal_to_bins(folder, path_to, totalbins, image_size_px, 
                                           image_size_um, camera_px_length, window, 
                                           deg, repetitions, mode, factor, meas_pow_bfp, 
                                           start_notch, end_notch, start_power, end_power, 
                                           start_spr, max_spr, lower_londa, upper_londa, 
                                           exp_time_slow, w0, use_calibration_flag, 
                                           last_bin_is_bkg_flag, 
                                           plot_flag=plot_flag)
    else:
        print('\nSTEP 1 was not executed.')
        
    #####################################################################
    #####################################################################
    #####################################################################
        
    if do_step2:
        for f in list_of_folders_slow:
            if single_NP_flag:    
                if not re.search('_NP_%03d' % NP_int, f):
                    continue
                else:
                    print('\nTarget NP found!')
                
            folder = os.path.join(parent_folder,f)
        
            print('\n === NP folder:',f)
        
            print('\n------------------------ STEP 2')

            step2.calculate_quotient(folder, path_to, totalbins, lower_londa, 
                                     upper_londa, Tzero_guess, beta_guess, 
                                     last_bin_is_bkg_flag, find_beta, find_Tzero)
                
    else:
        print('\nSTEP 2 was not executed.')
        
    #####################################################################
    #####################################################################
    #####################################################################
    
    if do_step3:
        for f in list_of_folders_slow:
            if single_NP_flag:    
                if not re.search('_NP_%03d' % NP_int, f):
                    continue
                else:
                    print('\nTarget NP found!')
                
            folder = os.path.join(parent_folder,f)
        
            print('\n === NP folder:',f)
            
            print('\n------------------------ STEP 3')
        
            step3.calculate_temp(folder, path_to, totalbins, alpha, R2th)
    else:
        print('\nSTEP 3 was not executed.')
    
    #####################################################################
    #####################################################################
    #####################################################################
    
    if do_step4 and not single_NP_flag:
    
        print('\n------------------------ STEP 4')
    
        step4.gather_data(path_to, R2th, totalbins, monitor_flag)
        
        step4.statistics(path_to, R2th, totalbins, radius, sigma_abs, find_beta, find_Tzero)
    else:
        print('\nSTEP 4 was not executed.')
    
    print('\nProcess done.')
    
    return

if __name__ == '__main__':
    
    beta_guess = 63.4 # beta bb for 80 AuNPs
#    beta_guess = 62.547 # beta bb for a single 80 AuNPs, the one scanned at room temp in the heating stage temp
#    beta_guess = 62.00 # beta bb for a single 80 AuNPs, averaged across all temp in the heating stage temp
    
#    T_lab_guess = 22
#    meas_pow_bfp = 0.646 # in mW
#    NP_folder = 'AuNP_SS_80/20191015_power/20191015-163050_Luminescence_10x10_power_0646uW'

#    T_lab_guess = 22
#    meas_pow_bfp = 0.42 # in mW
#    NP_folder = 'AuNP_SS_153/20191121_marienfeld/20191121-185854_Luminescence_10x10'
#
#    T_lab_guess = 22
#    meas_pow_bfp = 1.66 # in mW
#    NP_folder = 'AuNPs_SS_50/20191219-182356_Luminescence_10x10'
#
#    T_lab_guess = 22
#    meas_pow_bfp = 0.77 # in mW   
#    NP_folder = 'AuNPz_NPC_60magicas/set_de_mediciones_automaticas/20190902_repetitividad/20190902-164955_Luminescence 10x4_for_new_analysis'
#
#    T_lab_guess = 22
#    meas_pow_bfp = 0.671 # in mW
#    NP_folder = 'AuNP_SS_80/20191017_sapphire/20191017-205438_Luminescence_10x10'
    
#    T_lab_guess = 22
#    meas_pow_bfp = 0.652 # in mWu
#    NP_folder = 'AuNP_SS_80/20200312_grafeno/20200312-144420_Luminescence_10x10/'
#    NP_folder = 'AuNP_SS_80/20200312_grafeno/20200312-194257_Luminescence_10x6/'


#    Data analysis for different power
    T_lab_guess = 22
    meas_pow_bfp = 0.333 # in mW
    NP_folder = 'AuNP_SS_80/20191015_power/20191015-202833_Luminescence_10x10_power_0333uW'

#    T_lab_guess = 22
#    meas_pow_bfp = 1.013 # in mW
#    NP_folder = 'AuNP_SS_80/20191015_power/20191016-065230_Luminescence_10x10_power_1013uW'

    find_beta = True
    find_Tzero = False
#    
#    T_lab_guess = 23
#    NP_folder = 'AuNP_SS_80/20191123_temperatura/20191203-175127_Luminescence_10x4_22'
#    meas_pow_bfp = 0.630 # in mW
#
#    T_lab_guess = 30
#    NP_folder = 'AuNP_SS_80/20191123_temperatura/20191203-191827_Luminescence_11x1_30'
#    meas_pow_bfp = 0.630 # in mW
#
#    T_lab_guess = 40
#    NP_folder = 'AuNP_SS_80/20191123_temperatura/20191203-202747_Luminescence_11x1_40'
#    meas_pow_bfp = 0.630 # in mW
#
#    T_lab_guess = 50
#    NP_folder = 'AuNP_SS_80/20191123_temperatura/20191204-114857_Luminescence_11x1_50'
#    meas_pow_bfp = 0.622*1.01 # in mW, correction factor due to power monitor observation -> 0.628
###
#    T_lab_guess = 60
#    NP_folder = 'AuNP_SS_80/20191123_temperatura/20191204-125150_Luminescence_11x1_60'
#    meas_pow_bfp = 0.622*1.04 # in mW, correction factor due to power monitor observation -> 0.649
#    
#    T_lab_guess = 70
#    NP_folder = 'AuNP_SS_80/20191123_temperatura/20191204-155201_Luminescence_11x1_70'
#    meas_pow_bfp = 0.622*1.06 # in mW, correction factor due to power monitor observation -> 0.659

#    T_lab_guess = 25
#    NP_folder = 'AuNP_SS_80/20191123_temperatura/20191204-180125_Luminescence_11x1_25'
#    meas_pow_bfp = 0.670 # in mW
 
#    find_beta = False
#    find_Tzero = True
    
    print('\nProcessing %s' % NP_folder)
    print('\nPhotothermal coefficient (initial guess): %.1f K µm2/mW' % beta_guess) 
    print('\nLab. temperature (initial guess): %.1f °C' % T_lab_guess) 
    print('\nPower at BFP: %.3f mW' % meas_pow_bfp) 
    run_analysis(T_lab_guess, beta_guess, find_beta, find_Tzero, meas_pow_bfp, NP_folder)
        
        
        
        