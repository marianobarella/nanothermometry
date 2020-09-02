# -*- coding: utf-8 -*-
"""
Analysis of temperature increase of single AuNPs

Mariano Barella

21 aug 2018

CIBION, Buenos Aires, Argentina
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from functions_for_photoluminiscence import manage_save_directory

try:
    plt.style.use('for_confocal.mplstyle')
except:
    print('Pre-defined matplotlib style was not loaded.')

plt.ioff()
plt.close('all')

def calculate_temp(folder, path_to, totalbins, alpha, R2th):
    
    # Load folders
    NP =  folder.split('Spectrum_')[-1]
    
    save_folder = os.path.join(path_to, NP)
    
    common_path = os.path.join(path_to,'common_plots')

    folder = os.path.join(save_folder, 'matrix')
    
    bin_folder = os.path.join(save_folder, 'pl_in_bins')
    
    beta_matrix_file = os.path.join(folder, 'beta_matrix_%s.dat' % NP)
    beta_matrix = np.loadtxt(beta_matrix_file, delimiter=',')
    
    err_beta_matrix_file = os.path.join(folder, 'err_beta_matrix_%s.dat' % NP)
    err_beta_matrix = np.loadtxt(err_beta_matrix_file, delimiter=',')
    
    Tzero_matrix_file = os.path.join(folder, 'Tzero_matrix_%s.dat' % NP)
    Tzero_matrix = np.loadtxt(Tzero_matrix_file, delimiter=',')
    
    err_Tzero_matrix_file = os.path.join(folder, 'err_Tzero_matrix_%s.dat' % NP)
    err_Tzero_matrix = np.loadtxt(err_Tzero_matrix_file, delimiter=',')
    
    T_matrix_file = os.path.join(folder, 'Temp_matrix_%s.dat' % NP)
    T_matrix = np.loadtxt(T_matrix_file, delimiter=',')
    
    err_T_matrix_file = os.path.join(folder, 'err_T_matrix_%s.dat' % NP)
    err_T_matrix = np.loadtxt(err_T_matrix_file, delimiter=',')
    
    R2_matrix_file = os.path.join(folder, 'R2_matrix_%s.dat' % NP)
    R2_matrix = np.loadtxt(R2_matrix_file, delimiter=',')
    
    p_value_matrix_file = os.path.join(folder, 'p_value_matrix_%s.dat' % NP)
    p_value_matrix = np.loadtxt(p_value_matrix_file, delimiter=',')
    
    list_of_files = os.listdir(bin_folder)
    list_of_files.sort()
    corrected_mean_irrad_file = os.path.join(bin_folder,'bin_irradiance_%s.dat' % NP)
    data_irrad = np.loadtxt(corrected_mean_irrad_file)
    mean_irrad = data_irrad[:, 0]
    
    print('\n-- NP ', NP)

    # ALLOCATION
    T_good_matrix = np.zeros([totalbins,totalbins])
    err_T_good_matrix = np.zeros([totalbins,totalbins])
    Tzero_good_matrix = np.zeros([totalbins,totalbins])
    err_Tzero_good_matrix = np.zeros([totalbins,totalbins])
    beta_good_matrix = np.zeros([totalbins,totalbins])
    err_beta_good_matrix = np.zeros([totalbins,totalbins])
    irrad_good = np.array([])
    T_avg = np.array([])
    T_err = np.array([])
    
    # APPLY CRITERIA TO ELIMINATE UNDESIRED/OUTLIERS/BAD-FITTED POINTS
    good_p_value = p_value_matrix > alpha
    good_r2 = R2_matrix > R2th
    good = good_r2 & good_p_value
    
    T_good_matrix[good] = T_matrix[good]
    err_T_good_matrix[good] = err_T_matrix[good]
    
    Tzero_good_matrix[good] = Tzero_matrix[good]
    err_Tzero_good_matrix[good] = err_Tzero_matrix[good]
    
    beta_good_matrix[good] = beta_matrix[good]
    err_beta_good_matrix[good] = err_beta_matrix[good]
    
    ################# PLOT GOOD TEMP
    ################# PLOT GOOD TEMP
    ################# PLOT GOOD TEMP
    
    plt.figure()
    plt.imshow(T_good_matrix, interpolation='none', cmap='plasma')
    ax = plt.gca()
    for i in range(totalbins):
        for j in range(totalbins):
            ax.text(j, i, '%.0f' % T_good_matrix[i, j],
                    ha="center", va="center", color=[0,0,0])
    ax.set_xticks(range(totalbins))
    ax.set_yticks(range(totalbins))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    cbar = plt.colorbar()
    cbar.ax.set_title(u'Temp. (K)', fontsize=13)
    figure_name = os.path.join(folder,'temp_good_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'temp_matrix')
    figure_name = os.path.join(aux_folder,'temp_good_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()
    
    plt.figure()
    plt.imshow(Tzero_good_matrix, interpolation='none', cmap='plasma')
    ax = plt.gca()
    for i in range(totalbins):
        for j in range(totalbins):
            ax.text(j, i, '%.0f' % Tzero_good_matrix[i, j],
                    ha="center", va="center", color=[0,0,0])
    ax.set_xticks(range(totalbins))
    ax.set_yticks(range(totalbins))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    cbar = plt.colorbar()
    cbar.ax.set_title(u'T$_{0}$ (K)', fontsize=13)
    figure_name = os.path.join(folder,'Tzero_good_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'Tzero_matrix')
    figure_name = os.path.join(aux_folder,'Tzero_good_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()
    
    plt.figure()
    plt.imshow(beta_good_matrix, interpolation='none', cmap='plasma')
    ax = plt.gca()
    for i in range(totalbins):
        for j in range(totalbins):
            ax.text(j, i, '%.0f' % beta_good_matrix[i, j],
                    ha="center", va="center", color=[0,0,0])
    ax.set_xticks(range(totalbins))
    ax.set_yticks(range(totalbins))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    cbar = plt.colorbar()
    cbar.ax.set_title('ß (K µm$^{2}$/mW)', fontsize=13)
    figure_name = os.path.join(folder,'beta_good_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'beta_matrix')
    figure_name = os.path.join(aux_folder,'beta_good_matrix_%s.png' % NP)
    plt.savefig(figure_name)
    plt.close()
    
    ################# PLOT GOOD TZERO
    ################# PLOT GOOD TZERO
    ################# PLOT GOOD TZERO
    
    # calc Tzero
    Tzero_good_array = Tzero_good_matrix.reshape(1,totalbins**2)[0]
    Tzero_good_array = Tzero_good_array[Tzero_good_array != 0]
    
    err_Tzero_good_array = err_Tzero_good_matrix.reshape(1,totalbins**2)[0]
    err_Tzero_good_array = err_Tzero_good_array[err_Tzero_good_array != 0]
    
    mu_Tzero = np.mean(Tzero_good_array)
    std_Tzero = np.std(Tzero_good_array, ddof = 1)
    borders_Tzero = [mu_Tzero + std_Tzero, mu_Tzero - std_Tzero]
#    err_mu_Tzero = np.sqrt(np.sum(err_Tzero_good_array**2))/len(Tzero_good_array)
    err_mu_Tzero = std_Tzero/np.sqrt(len(Tzero_good_array))
    
    plt.figure()
    nbins = 20
#    range_tuple = [np.min(beta_good_array), np.max(beta_good_array)]
    range_tuple = [0,100]
    plt.hist(Tzero_good_array, bins=nbins, range=range_tuple, rwidth = 1, \
             align='mid', color='C0', alpha = 1, edgecolor='k', normed = False)
#    plt.legend(loc=1, handlelength=0)
    ax = plt.gca()
    ax.axvline(mu_Tzero, ymin = 0, ymax = 1, color='k', linestyle='--', linewidth=2)
    ax.fill_between([borders_Tzero[0], borders_Tzero[1]], [0,0], [100,100],  facecolor='k', alpha=0.25)
#    ax.get_yaxis().set_ticks([])
    ax.set_ylim([0, 10])
    ax.set_xlim([1, 99])
    ax.set_xlabel(u'T$_{0}$ (K)')
    ax.set_ylabel('Entries')
    aux_folder = manage_save_directory(folder,'Tzero_matrix')
    figure_name = os.path.join(aux_folder,'Tzero_good_hist_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'Tzero_matrix')
    figure_name = os.path.join(aux_folder,'Tzero_good_hist_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300, bbox_layout = 'tight')
    plt.close()
    
    ################# PLOT GOOD BETA
    ################# PLOT GOOD BETA
    ################# PLOT GOOD BETA
    
    # calc beta
    beta_good_array = beta_good_matrix.reshape(1,totalbins**2)[0]
    beta_good_array = beta_good_array[beta_good_array != 0]
    
    err_beta_good_array = err_beta_good_matrix.reshape(1,totalbins**2)[0]
    err_beta_good_array = err_beta_good_array[err_beta_good_array != 0]

    data_points = len(beta_good_array)
    mu_beta = np.mean(beta_good_array)
    std_beta = np.std(beta_good_array, ddof = 1)
    borders_beta = [mu_beta + std_beta, mu_beta - std_beta]
#    err_mu_beta = np.sqrt(np.sum(err_beta_good_array**2))/len(beta_good_array)
    err_mu_beta = std_beta/np.sqrt(data_points)
    
#    print(mu_beta, std_beta, err_mu_beta)
    plt.figure()
    nbins = 20
#    range_tuple = [np.min(beta_good_array), np.max(beta_good_array)]
    range_tuple = [0,100]
    plt.hist(beta_good_array, bins=nbins, range=range_tuple, rwidth = 1, \
             align='mid', color='C0', alpha = 1, edgecolor='k', normed = False)
#    plt.legend(loc=1, handlelength=0)
    ax = plt.gca()
    ax.axvline(mu_beta, ymin = 0, ymax = 1, color='k', linestyle='--', linewidth=2)
    ax.fill_between([borders_beta[0], borders_beta[1]], [0,0], [100,100],  facecolor='k', alpha=0.25)
#    ax.get_yaxis().set_ticks([])
    ax.set_ylim([0, 10])
    ax.set_xlim([1, 99])
    ax.set_xlabel(u'Photothermal coefficient (K µm$^{2}$/mW)')
    ax.set_ylabel('Entries')
    aux_folder = manage_save_directory(folder,'beta_matrix')
    figure_name = os.path.join(aux_folder,'beta_good_hist_%s.png' % NP)
    plt.savefig(figure_name)
    aux_folder = manage_save_directory(common_path,'beta_matrix')
    figure_name = os.path.join(aux_folder,'beta_good_hist_%s.png' % NP)
    plt.savefig(figure_name, dpi = 300, bbox_layout = 'tight')
    plt.close()
    
    ########################### FIND SLOPE/BETA FITTING TEMPS
    ########################### FIND SLOPE/BETA FITTING TEMPS
    ########################### FIND SLOPE/BETA FITTING TEMPS
    
    # Do some minor statistics
    for i in range(totalbins):
        N = np.sum(good[i,:])
        if not N: 
            print('Bin %d has no temperature values that fulfill our criteria.' % i)
            continue
        print('Bin %d has %d temperature values.' % (i, N))
        # temp list
        index_T_good = np.where(T_good_matrix[i,:] != 0)[0]
        T_list = T_good_matrix[i,index_T_good]
        T_list = np.array(T_list)
        # error list
        err_T_list = err_T_good_matrix[i,index_T_good]
        err_T_list = np.array(err_T_list)
        # calculate mean and error
        aux_mean = np.sum(T_list)/N
        T_avg = np.append(T_avg, aux_mean)
        # append error of average of T_list (matrix's row)
        aux_err = np.sqrt(np.sum(err_T_list**2))/N
        T_err = np.append(T_err, aux_err)
        # append irradiance
        irrad_good = np.append(irrad_good, mean_irrad[i])
#        print('Temp. increase', T_list)
#        print('Errors', err_T_list)
#        print('Weights', list_T_weights)
       
   	# Fit slope for each NP in order to find information about medium dissipation
    x = irrad_good
    y = T_avg
    y_err = T_err
	
    ############### PLOT #####################
    ############### PLOT #####################
    ############### PLOT #####################
    x_beta = np.arange(0, 5, 0.1)
    y_beta = x_beta*mu_beta + mu_Tzero

    plt.figure()
    plt.errorbar(x, y, yerr = y_err, fmt = 'o', color = 'C0',
                 ms = 7, mfc = 'C0', ecolor = 'C0', lw = 1, capsize = 2.5, 
                 barsabove = False)
    plt.plot(x_beta, y_beta, 'k-')
    ax = plt.gca()
    ax.set_xlabel(u'Irradiance (mW/µm$^{2}$)', fontsize=20)
    ax.set_ylabel('Temperature (K)', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=18)
    ax.set_xlim([0, 2])
#    x_axis = list(np.arange(0,3.1,0.5))
#    ax.set_xticks(x_axis)
#    ax.set_xticklabels(x_axis)
    ax.set_ylim([290, 400])
#        ax.set_ylim([-10, 150])
    aux_folder = manage_save_directory(folder,'temp_vs_irrad')
    figure_name = os.path.join(aux_folder, 'temp_vs_irrad_R2th_%s_%s.png' % (str(R2th),NP))
    plt.savefig(figure_name, dpi=300)
    aux_folder = manage_save_directory(common_path,'temp_vs_irrad')
    figure_name = os.path.join(aux_folder, 'temp_vs_irrad_R2th_%s_%s.png' % (str(R2th),NP))
    plt.savefig(figure_name, dpi=300)
    figure_name = os.path.join(aux_folder, 'temp_vs_irrad_R2th_%s_%s.pdf' % (str(R2th),NP))
    plt.savefig(figure_name, dpi=300, format = 'pdf')
    plt.close()

    #################### SAVE DATA ############################
    #################### SAVE DATA ############################
    #################### SAVE DATA ############################

    save_to_array = np.array([x, y, y_err]).T
    aux_folder = manage_save_directory(folder,'temp_vs_irrad')
    temp_vs_irrad_file = os.path.join(aux_folder,'temp_vs_irrad_R2th_%s_%s.txt' % (str(R2th),NP))
    np.savetxt(temp_vs_irrad_file, save_to_array, delimiter=',', fmt='%.3e')
    
    aux_folder = manage_save_directory(folder,'beta_matrix')
    slope_file = os.path.join(aux_folder, 'beta_R2th_%s_%s.dat' % (str(R2th),NP))
    f = open(slope_file, 'w+')
    f.write('MEAN_BETA_(K*um2/mW) ERROR_OF_MEAN_BETA_(K*um2/mW) STD_BETA_(K*um2/mW) DATA_POINTS MEAN_TZERO_(K) ERROR_OF_MEAN_TZERO_(K) STD_TZERO_(K)\n')
    string_to_write = '%.3e %.3e %.3e %d %.3e %.3e %.3e\n' % (mu_beta, err_mu_beta, std_beta, data_points, mu_Tzero, err_mu_Tzero, std_Tzero)
    f.write(string_to_write)
    f.close()
    
    beta_good_array_file = os.path.join(aux_folder,'beta_good_array_%s.dat' % NP)
    np.savetxt(beta_good_array_file, beta_good_array, delimiter=',', fmt='%.3e')
    
    err_beta_good_array_file = os.path.join(aux_folder,'beta_good_err_array_%s.dat' % NP)
    np.savetxt(err_beta_good_array_file, err_beta_good_array, delimiter=',', fmt='%.3e')
    
    aux_folder = manage_save_directory(folder,'Tzero_matrix')
    
    Tzero_good_array_file = os.path.join(aux_folder,'Tzero_good_array_%s.dat' % NP)
    np.savetxt(Tzero_good_array_file, Tzero_good_array, delimiter=',', fmt='%.3e')
    
    err_Tzero_good_array_file = os.path.join(aux_folder,'Tzero_good_err_array_%s.dat' % NP)
    np.savetxt(err_Tzero_good_array_file, err_Tzero_good_array, delimiter=',', fmt='%.3e')
    
    return