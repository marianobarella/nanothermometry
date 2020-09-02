# -*- coding: utf-8 -*-
"""
Process of single AuNPs spectral confocal image 
acquired with PySpectrum at CIBION

Mariano Barella

16 aug 2019

based on "witec_data_photoluminiscense.py"

"""

import os
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import cm
import scipy
from functions_for_photoluminiscence import manage_save_directory,\
    lambda_to_energy, calc_r2, calc_chi_squared, closer, skip_neighbours, \
    temp_with_slope, quotient_with_slope, fit_quotient_for_beta, \
    fit_quotient_for_Tzero, bose, log_bose, fit_bose, fit_log_bose

try:
    plt.style.use('for_confocal.mplstyle')
except:
    print('Pre-defined matplotlib style was not loaded.')

plt.ioff()
plt.close('all')

def calculate_quotient(folder, path_to, totalbins, lower_londa, 
                       upper_londa, Tzero_guess, beta_guess, 
                       last_bin_is_bkg_flag, find_beta, find_Tzero):

    thermometry_using_scattering = True
        
    NP = folder.split('Spectrum_')[-1]
    
    common_path = os.path.join(path_to,'common_plots')
    
    save_folder = os.path.join(path_to, NP)
    
    bin_folder = os.path.join(save_folder, 'pl_in_bins')
    
    londa_file = os.path.join(bin_folder,'londa_%s.dat' % NP)
    londa = np.loadtxt(londa_file)
    lower_londa_index = closer(londa, lower_londa)
    upper_londa_index = closer(londa, upper_londa)
    energy = lambda_to_energy(londa)
    
    a = lower_londa_index
    b = upper_londa_index
        
    dof = b - a + 1 - 1 # número de datos (ver ajuste con energía) MENOS los parámetros del ajuste 
    frozen_distro = scipy.stats.chi2(dof)
    
    bin_specs_file = os.path.join(bin_folder,'all_bins_%s.dat' % NP)
    bin_specs = np.loadtxt(bin_specs_file)
    bin_specs = bin_specs.T
    
    corrected_mean_irrad_file = os.path.join(bin_folder,'bin_irradiance_%s.dat' % NP)
    irrad = np.loadtxt(corrected_mean_irrad_file)
    Irrad_array = irrad[:,0]
    err_Irrad_array = irrad[:,1]
    
    bin_distro_file = os.path.join(bin_folder,'bin_distro_%s.dat' % NP)
    bin_data = np.loadtxt(bin_distro_file)
    hist_bin = bin_data[:,2]
        
    R2_matrix = np.zeros([totalbins, totalbins])
    chi_squared_matrix = np.zeros([totalbins, totalbins])
    p_value_matrix = np.zeros([totalbins, totalbins])
    beta_matrix = np.zeros([totalbins, totalbins])
    err_beta_matrix = np.zeros([totalbins, totalbins])
    Tzero_matrix = np.zeros([totalbins, totalbins])
    err_Tzero_matrix = np.zeros([totalbins, totalbins])
    T_matrix = np.zeros([totalbins, totalbins])
    err_T_matrix = np.zeros([totalbins, totalbins])
        
    if last_bin_is_bkg_flag:
        print('Last bin is considered a bkg correction, quotients calculated from bin 1.')
        list_of_bins = list(range(1, totalbins))
    else:
        print('Last bin is NOT a bkg correction, quotients calculated from bin 0.')
        list_of_bins = list(range(0, totalbins))
     
    xdata_energy = energy[a:b]
    xdata_londa = londa[a:b]
    len_data = len(xdata_energy)
        
    print('Allocating, calculating and plotting quotients')
    # allocate arrays for experimental data
    bin_A_array = np.zeros((totalbins, totalbins))
    err_bin_A_array = np.zeros((totalbins, totalbins))
    ydata_all = np.zeros((totalbins, totalbins, len_data))
    
    for ref in list_of_bins:
        if hist_bin[ref] == 0:
            print('\nBin %d is empty.' % ref)
            print('Skipping bin as reference.')
            continue
        # define reference irradiance
        Irrad_ref = Irrad_array[ref]
        err_Irrad_ref = err_Irrad_array[ref]
        # define quotients of irradiance
        for i in list_of_bins:    
            A = Irrad_array[i]/Irrad_ref
            bin_A_array[i, ref] = A
            err_bin_A_array[i, ref] = np.sqrt( (err_Irrad_array[i]/Irrad_ref)**2 + \
                                               (err_Irrad_ref*Irrad_array[i]/(Irrad_ref**2))**2)
            # load experimental data
            quotient_exp = bin_specs[i][a:b]/bin_specs[ref][a:b]
            ydata_all[i,ref,:] = quotient_exp
            
    if thermometry_using_scattering:
        print('\n Warning: Finding beta (also) using scattering spectrum!')
        # load scattering spectrum
        filepath_scattering_80 = '/home/mariano/datos_mariano/posdoc/experimentos_PL_arg/20190717_espectros_scattering_PL/grilla3_AuNPs_80SS/figures/averaged_scattering_spectra_80nm_glass_water.dat'
        data = np.loadtxt(filepath_scattering_80)
        londa_sca = data[:,0]
        spec_sca = data[:,1]
        energy_sca = lambda_to_energy(londa_sca)
    
        # Interpolate
        # note that inversion of arrays is needed to interpolate correctly
        new_sca_spectrum = np.interp(xdata_energy, energy_sca[::-1], spec_sca[::-1]) 
        # CHECK
        plt.figure()
        plt.plot(energy_sca, spec_sca,'C0-')
        plt.plot(xdata_energy, new_sca_spectrum,'k-', linewidth=4, alpha=0.4)
        ax = plt.gca()
        ax.set_xlabel('Energy (eV)', fontsize=20)
        ax.set_ylabel('Normalized scattering (a.u.)', fontsize=20)
        ax.set_xlim([1.71, 2.59])
        ax.set_ylim([0, 1.025])
        plt.tick_params(axis='both', which='major', labelsize=18)
        aux_folder = manage_save_directory(common_path,'antistokes_sca_and_pl')
        figure_name = os.path.join(aux_folder,'sca_interpolation_%s.png' % (NP))
        plt.savefig(figure_name, dpi = 300)
        figure_name = os.path.join(aux_folder,'sca_interpolation_%s.pdf' % (NP))
        plt.savefig(figure_name, dpi = 300, format='pdf')
        plt.close() 
        
        x_long = np.arange(2.37, 2.47, 0.01)    
        rango = list(range(0, totalbins + 1))
        norm = plt.Normalize()
        colormap = plt.cm.plasma(norm(rango))  
        
        temp_array_with_sca = np.zeros((totalbins))
        temp_array_with_sca[0] = Tzero_guess
        err_temp_array_with_sca = np.zeros((totalbins))
        err_temp_array_with_sca[0] = 1
        
        plt.figure()
#        plt.plot(xdata_energy, log_bose(xdata_energy, 1, 300),'r-')
#        plt.plot(xdata_energy, new_sca_spectrum,'.-')
        for i in list_of_bins:
            if hist_bin[i] == 0:
                print('\nBin %d is empty.' % i)
                print('Skipping bin.')
                continue
            # find temp or beta
            cte = 1
            y = cte*bin_specs[i][a:b]/new_sca_spectrum
            log_y = np.log10(y)
            init_params = [1, Tzero_guess]
            best_as, err = fit_log_bose(init_params, xdata_energy, log_y)
            temp_array_with_sca[i] = best_as[1]
            err_temp_array_with_sca[i] = np.sqrt(err[1,1])
            fitted_log_y = 10**(log_bose(x_long, best_as[0], best_as[1]))/1000
            plt.plot(x_long, fitted_log_y, 'k--')
            color_iter = colormap[i]
            plt.plot(xdata_energy, 10**(log_y)/1000, '-', linewidth = 2.0,
                         color = color_iter,
                         label='%d, T=%.0f K' % (i, best_as[1]))
            
        plt.grid(True)
        plt.legend(loc='upper right', prop={'size':12})
        ax = plt.gca()
#        ax.set_xlim([energy[b],energy[a]])
#        ax.set_xlim([0.99*energy[b],1.01*energy[a]])
        ax.set_xlim(2.37, 2.445)
#        ax.axvline(energy[a], ymin=0, ymax=1, linestyle='--', color='k')
#        ax.axvline(energy[b], ymin=0, ymax=1, linestyle='--', color='k')
        ax.set_ylim([0.1, 20])
        ax.set_xlabel('Energy (eV)', fontsize=20)
        ax.set_ylabel('AS/scattering (a.u.)', fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=18)
        aux_folder = manage_save_directory(common_path,'antistokes_sca_and_pl')
        figure_name = os.path.join(aux_folder,'pl_and_sca_comparison_%s.png' % (NP))
        plt.savefig(figure_name, dpi = 300)
        figure_name = os.path.join(aux_folder,'pl_and_sca_comparison_%s.pdf' % (NP))
        plt.savefig(figure_name, dpi = 300, format='pdf')
        plt.close()

        # plot temp vs irrad
        plt.figure()
        plt.errorbar(Irrad_array, temp_array_with_sca, 
                     xerr = err_Irrad_array,
                     yerr = err_temp_array_with_sca,
                     fmt = 'o', color = 'C0',
                     ms = 7, mfc = 'C0', ecolor = 'C0', lw = 1, capsize = 2.5, 
                     barsabove = False)
        plt.grid(True)
        plt.legend(loc='upper right', prop={'size':12})
        ax = plt.gca()
        ax.set_xlim([-0.05, 1.75])
        ax.set_ylim([285, 1125])
        ax.set_xlabel(u'Irradiance (mW/µm$^{2}$)', fontsize=20)
        ax.set_ylabel('Temperature (K)', fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=18)
        aux_folder = manage_save_directory(common_path,'antistokes_sca_and_pl')
        figure_name = os.path.join(aux_folder,'temp_vs_irrad_%s.png' % (NP))
        plt.savefig(figure_name, dpi = 300)
        figure_name = os.path.join(aux_folder,'temp_vs_irrad_%s.pdf' % (NP))
        plt.savefig(figure_name, dpi = 300, format='pdf')
        plt.close()   
            
    ##################### ANTI-STOKES vs LAMBDA ######################
    ##################### ANTI-STOKES vs LAMBDA ######################
    ##################### ANTI-STOKES vs LAMBDA ######################
    
    for ref in list_of_bins:
        plt.figure()
        for i in list_of_bins:    
            plt.plot(xdata_londa, ydata_all[i,ref,:], '-', label='Bin %d/%d - A %.2f' % (i, ref, bin_A_array[i, ref]))
        plt.grid(True)
        plt.legend(loc='best')
        ax = plt.gca()
        ax.set_xlabel(r'Wavelength (nm)')
        ax.set_ylabel('$Q^{AS}$')
        ax.set_xlim([0.99*londa[a],1.01*londa[b]])
        ax.axvline(londa[a], ymin=0, ymax=1, linestyle='--', color='k')
        ax.axvline(londa[b], ymin=0, ymax=1, linestyle='--', color='k')
        aux_folder = manage_save_directory(save_folder,'antistokes_quotient')
        figure_name = os.path.join(aux_folder,'quotient_vs_lambda_ref_%02d_%s.png' % (ref, NP))
        plt.savefig(figure_name)
        plt.close()

    ##################### ANTI-STOKES FIT ######################
    ##################### ANTI-STOKES FIT ######################
    ##################### ANTI-STOKES FIT ######################                
    x = energy[a:b]
    x_long = np.arange(2.37, 2.47, 0.01)    
    rango = list(range(0, totalbins + 1))
    norm = plt.Normalize()
    colormap = plt.cm.plasma(norm(rango))  
    
    for ref in list_of_bins:
        print('\nFitting quotient Q using bin %d as reference...' % ref)
        plt.figure()
        list_of_bins_to_fit = list(range(ref + 1 + skip_neighbours, totalbins))
        for i in list_of_bins_to_fit:
            # m = ref
            Irrad_n = Irrad_array[i]
            Irrad_m = Irrad_array[ref]
            err_Irrad_n = err_Irrad_array[i]
            err_Irrad_m = err_Irrad_array[ref]
            y = ydata_all[i, ref, :]
            if find_beta and not find_Tzero:
                print('\nTzero is known. Finding beta...')
                Tzero = Tzero_guess
                init_params = beta_guess # initial guess for beta
                # Get the fitting parameters for the best quotient of photoluminiscence emission
                best_as, err = fit_quotient_for_beta(init_params, x, Tzero, Irrad_n, Irrad_m, y)
                # retrieve best fitted parameters
                beta = best_as[0]
                Tref = temp_with_slope(Tzero, beta, Irrad_m)
                Tbin = temp_with_slope(Tzero, beta, Irrad_n)
                # calculate the errors and goodes of the fit
                err_beta = np.sqrt(err[0,0])
                err_Tzero = 0.5
                err_Tref = np.sqrt(err_Tzero**2 + (Irrad_m*err_beta)**2 + (err_Irrad_m*beta)**2)
                err_Tbin = np.sqrt(err_Tzero**2 + (Irrad_n*err_beta)**2 + (err_Irrad_n*beta)**2)
                yfit = quotient_with_slope(x, Tzero, beta, Irrad_n, Irrad_m)
                r2_coef_pearson = calc_r2(y, yfit)
                chi_squared_pearson = calc_chi_squared(y, yfit)
                # Plotting
                yfit_to_plot = quotient_with_slope(x_long, Tzero, beta, Irrad_n, Irrad_m)
                plt.plot(x_long, yfit_to_plot, '--k', alpha = 0.8)
                color_iter = colormap[i]
                plt.plot(x, y, '-', linewidth = 2.0,
                         color = color_iter,
                         label='$Q_{%d/%d}$' % (i, ref))
    #                     label='$Q^{AS}_{%d/%d}$' % (i, ref))
    #                     label='Group %d/%d' % (i, ref))
    #                    label='Bin %d/%d - A %.2f' % (i, ref, A))
                # Asign matrix elements
                beta_matrix[ref,i] = beta
                err_beta_matrix[ref,i] = err_beta
                Tzero_matrix[ref,i] = Tzero
                err_Tzero_matrix[ref,i] = err_Tzero
                T_matrix[ref,i] = Tref
                T_matrix[i,ref] = Tbin
                err_T_matrix[ref,i] = err_Tref
                err_T_matrix[i,ref] = err_Tbin
                R2_matrix[ref,i] = r2_coef_pearson
                R2_matrix[i,ref] = r2_coef_pearson
                chi_squared_matrix[ref,i] = chi_squared_pearson
                chi_squared_matrix[i,ref] = chi_squared_pearson
                p_value_matrix[ref,i] = 1 - frozen_distro.cdf(chi_squared_pearson)
                p_value_matrix[i,ref] = 1 - frozen_distro.cdf(chi_squared_pearson)
                print('---------- Bin', i, \
                      '\nR-sq: %.3f' % r2_coef_pearson, \
                      '\nTzero: %.1f' % Tzero, 'err Tzero: %.1f' % err_Tzero, \
                      '\nbeta: %.1f' % beta, 'error beta: %.2f' % err_beta, \
                      '\nT_ref: %.1f' % Tref, 'error T_ref: %.1f' % err_Tref, \
                      '\nT_bin: %.1f' % Tbin, 'error T_bin: %.1f' % err_Tbin)
                
            elif not find_beta and find_Tzero:
                print('\nBeta is known. Finding Tzero...')
                beta = beta_guess
                init_params = Tzero_guess # initial guess for Tzero
                # Get the fitting parameters for the best quotient of photoluminiscence emission
                best_as, err = fit_quotient_for_Tzero(init_params, x, beta, Irrad_n, Irrad_m, y)
                # retrieve best fitted parameters
                Tzero = best_as[0]
                Tref = temp_with_slope(Tzero, beta, Irrad_m)
                Tbin = temp_with_slope(Tzero, beta, Irrad_n)
                # calculate the errors and goodes of the fit
                err_beta = 0.05
                err_Tzero = np.sqrt(err[0,0])
                err_Tref = np.sqrt(err_Tzero**2 + (Irrad_m*err_beta)**2 + (err_Irrad_m*beta)**2)
                err_Tbin = np.sqrt(err_Tzero**2 + (Irrad_n*err_beta)**2 + (err_Irrad_n*beta)**2)
                yfit = quotient_with_slope(x, Tzero, beta, Irrad_n, Irrad_m)
                r2_coef_pearson = calc_r2(y, yfit)
                chi_squared_pearson = calc_chi_squared(y, yfit)
                # Plotting
                yfit_to_plot = quotient_with_slope(x_long, Tzero, beta, Irrad_n, Irrad_m)
                plt.plot(x_long, yfit_to_plot, '--k', alpha = 0.8)
                color_iter = colormap[i]
                plt.plot(x, y, '-', linewidth = 2.0,
                         color = color_iter,
                         label='$Q_{%d/%d}$' % (i, ref))
    #                     label='$Q^{AS}_{%d/%d}$' % (i, ref))
    #                     label='Group %d/%d' % (i, ref))
    #                    label='Bin %d/%d - A %.2f' % (i, ref, A))
                # Asign matrix elements
                beta_matrix[ref,i] = beta
                err_beta_matrix[ref,i] = err_beta
                Tzero_matrix[ref,i] = Tzero
                err_Tzero_matrix[ref,i] = err_Tzero
                T_matrix[ref,i] = Tref
                T_matrix[i,ref] = Tbin
                err_T_matrix[ref,i] = err_Tref
                err_T_matrix[i,ref] = err_Tbin
                R2_matrix[ref,i] = r2_coef_pearson
                R2_matrix[i,ref] = r2_coef_pearson
                chi_squared_matrix[ref,i] = chi_squared_pearson
                chi_squared_matrix[i,ref] = chi_squared_pearson
                p_value_matrix[ref,i] = 1 - frozen_distro.cdf(chi_squared_pearson)
                p_value_matrix[i,ref] = 1 - frozen_distro.cdf(chi_squared_pearson)
                print('---------- Bin', i, \
                      '\nR-sq: %.3f' % r2_coef_pearson, \
                      '\nTzero: %.1f' % Tzero, 'err Tzero: %.1f' % err_Tzero, \
                      '\nbeta: %.1f' % beta, 'error beta: %.2f' % err_beta, \
                      '\nT_ref: %.1f' % Tref, 'error T_ref: %.1f' % err_Tref, \
                      '\nT_bin: %.1f' % Tbin, 'error T_bin: %.1f' % err_Tbin)
            else:
                raise ValueError('Values of find_beta and find_Tzero must be complementary.')
                
        plt.grid(True)
        plt.legend(loc='upper left', prop={'size':12})
        ax = plt.gca()
        ax.set_xlim([energy[b],energy[a]])
#        ax.set_xlim([0.99*energy[b],1.01*energy[a]])
        ax.set_xlim(2.37, 2.445)
#        ax.axvline(energy[a], ymin=0, ymax=1, linestyle='--', color='k')
#        ax.axvline(energy[b], ymin=0, ymax=1, linestyle='--', color='k')
        ax.set_ylim([1, 24])
        ax.set_xlabel('Energy (eV)', fontsize=20)
        ax.set_ylabel('$Q^{AS}$', fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=18)
        aux_folder = manage_save_directory(save_folder,'antistokes_quotient')
        figure_name = os.path.join(aux_folder,'quotient_vs_energy_ref_%02d_%s.png' % (ref, NP))
        plt.savefig(figure_name, dpi = 300)
        figure_name = os.path.join(aux_folder,'quotient_vs_energy_ref_%02d_%s.pdf' % (ref, NP))
        plt.savefig(figure_name, dpi = 300, format='pdf')
        plt.close()
    
    ########################### PLOT Temp, R2, p-value matrix
    ########################### PLOT Temp, R2, p-value matrix
    
    plt.figure()
    plt.imshow(beta_matrix, interpolation='none', cmap='plasma')
    ax = plt.gca()
    for i in range(totalbins):
        for j in range(totalbins):
            ax.text(j, i, '%.0f' % beta_matrix[i, j],
                    ha="center", va="center", color=[0,0,0])
    ax.set_xticks(range(totalbins))
    ax.set_yticks(range(totalbins))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    cbar = plt.colorbar()
    cbar.ax.set_title(u'Temp. (K)', fontsize=13)
    figure_name = os.path.join(folder,'beta_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'beta_matrix')
    figure_name = os.path.join(aux_folder,'beta_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()
    
    plt.figure()
    plt.imshow(Tzero_matrix, interpolation='none', cmap='plasma')
    ax = plt.gca()
    for i in range(totalbins):
        for j in range(totalbins):
            ax.text(j, i, '%.0f' % Tzero_matrix[i, j],
                    ha="center", va="center", color=[0,0,0])
    ax.set_xticks(range(totalbins))
    ax.set_yticks(range(totalbins))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    cbar = plt.colorbar()
    cbar.ax.set_title(u'T$_{0}$ (K)', fontsize=13)
    figure_name = os.path.join(folder,'Tzero_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'Tzero_matrix')
    figure_name = os.path.join(aux_folder,'Tzero_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()
    
    plt.figure()
    plt.imshow(T_matrix, interpolation='none', cmap='plasma')
    ax = plt.gca()
    for i in range(totalbins):
        for j in range(totalbins):
            ax.text(j, i, '%.0f' % T_matrix[i, j],
                    ha="center", va="center", color=[0,0,0])
    ax.set_xticks(range(totalbins))
    ax.set_yticks(range(totalbins))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    cbar = plt.colorbar()
    cbar.ax.set_title(u'Temp. (K)', fontsize=13)
    figure_name = os.path.join(folder,'temp_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'temp_matrix')
    figure_name = os.path.join(aux_folder,'temp_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()
    
    plt.figure()
    plt.imshow(R2_matrix, interpolation='none', cmap='viridis')
    ax = plt.gca()
    for i in range(totalbins):
        for j in range(totalbins):
            ax.text(j, i, '%.2f' % R2_matrix[i, j],
                    ha="center", va="center", color=[0,0,0])
    ax.set_xticks(range(totalbins))
    ax.set_yticks(range(totalbins))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    figure_name = os.path.join(folder,'R2_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'R2_matrix')
    figure_name = os.path.join(aux_folder,'R2_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()
    
    plt.figure()
    plt.imshow(p_value_matrix, interpolation='none', cmap='viridis')
    ax = plt.gca()
    for i in range(totalbins):
        for j in range(totalbins):
            ax.text(j, i, '%.2f' % p_value_matrix[i, j],
                    ha="center", va="center", color=[0,0,0])
    ax.set_xticks(range(totalbins))
    ax.set_yticks(range(totalbins))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    figure_name = os.path.join(folder,'p_value_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'p_value_matrix')
    figure_name = os.path.join(aux_folder,'p_value_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()
    
    #################### SAVE DATA ############################
    #################### SAVE DATA ############################
    #################### SAVE DATA ############################
    
    # Writing/creating files  
    aux_folder = manage_save_directory(save_folder,'matrix')
    
    beta_matrix_file = os.path.join(aux_folder,'beta_matrix_%s.dat' % NP)
    err_beta_matrix_file = os.path.join(aux_folder,'err_beta_matrix_%s.dat' % NP)
    np.savetxt(beta_matrix_file, beta_matrix, delimiter=',', fmt='%.3e')
    np.savetxt(err_beta_matrix_file, err_beta_matrix, delimiter=',', fmt='%.3e')
    
    Tzero_matrix_file = os.path.join(aux_folder,'Tzero_matrix_%s.dat' % NP)
    err_Tzero_matrix_file = os.path.join(aux_folder,'err_Tzero_matrix_%s.dat' % NP)
    np.savetxt(Tzero_matrix_file, Tzero_matrix, delimiter=',', fmt='%.3e')
    np.savetxt(err_Tzero_matrix_file, err_Tzero_matrix, delimiter=',', fmt='%.3e')
    
    T_matrix_file = os.path.join(aux_folder,'Temp_matrix_%s.dat' % NP)
    err_T_matrix_file = os.path.join(aux_folder,'err_T_matrix_%s.dat' % NP)
    np.savetxt(T_matrix_file, T_matrix, delimiter=',', fmt='%.3e')
    np.savetxt(err_T_matrix_file, err_T_matrix, delimiter=',', fmt='%.3e')

    R2_matrix_file = os.path.join(aux_folder,'R2_matrix_%s.dat' % NP)
    np.savetxt(R2_matrix_file, R2_matrix, delimiter=',', fmt='%.3e')
    
    p_value_matrix_file = os.path.join(aux_folder,'p_value_matrix_%s.dat' % NP)
    np.savetxt(p_value_matrix_file, p_value_matrix, delimiter=',', fmt='%.3e')
            
    irradiance_file = os.path.join(aux_folder,'irradiance_matrix_%s.dat' % NP)  
    np.savetxt(irradiance_file, Irrad_array, fmt='%.3e') 
    
    return