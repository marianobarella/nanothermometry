# -*- coding: utf-8 -*-
"""
Analysis of temperature increase of single AuNPs

Mariano Barella

21 aug 2019

CIBION, Buenos Aires, Argentina
"""

import os
import re
import numpy as np
import pickle
import matplotlib.pyplot as plt
from functions_for_photoluminiscence import manage_save_directory

try:
    plt.style.use('for_confocal.mplstyle')
except:
    print('Pre-defined matplotlib style was not loaded.')

plt.ioff()
plt.close('all')

def gather_data(path_from, R2th, totalbins, monitor):

    #################### SLOPE STATISTICS ############################
    #################### SLOPE STATISTICS ############################
    #################### SLOPE STATISTICS ############################
    
    list_of_folders = os.listdir(path_from)
    list_of_folders = [f for f in list_of_folders if os.path.isdir(os.path.join(path_from,f))]
    list_of_folders = [f for f in list_of_folders if re.search('NP',f)]
    list_of_folders.sort()
    L = len(list_of_folders)
    
    list_of_mean_beta = np.zeros((L))
    list_of_mean_beta_err = np.zeros((L))
    list_of_std_beta = np.zeros((L))
    list_of_data_points = np.zeros((L))
    list_of_mean_Tzero = np.zeros((L))
    list_of_mean_Tzero_err = np.zeros((L))
    list_of_std_Tzero = np.zeros((L))
    list_of_londa_max = np.zeros((L))
    list_of_width_pl = np.zeros((L))
    list_of_r_sq_spr = np.zeros((L))
    list_of_col = np.zeros((L))
    list_of_nps = np.zeros((L))
    list_of_avg_monitor = np.zeros((L))
    list_of_std_monitor = np.zeros((L))
    
    list_of_beta_good_array = []
    list_of_beta_good_dict = {}
    list_of_beta_good_err_dict = {}
    
    list_of_Tzero_good_array = []
    list_of_Tzero_good_dict = {}
    list_of_Tzero_good_err_dict = {}
    
    for i in range(L):
        NP = list_of_folders[i]
    
        list_of_col[i] = int(NP.split('_')[1])
        list_of_nps[i] = int(NP.split('_')[3])
    
        ############################################ import beta
        beta_folder = os.path.join(path_from, NP, 'matrix', 'beta_matrix')
        filename = 'beta_R2th_%s_%s.dat' % (str(R2th), NP)
        beta_filepath = os.path.join(beta_folder, filename)
        data = np.loadtxt(beta_filepath, skiprows=1)
        list_of_mean_beta[i] = data[0]
        list_of_mean_beta_err[i] = data[1]
        list_of_std_beta[i] = data[2]
        list_of_data_points[i] = data[3]
        list_of_mean_Tzero[i] = data[4]
        list_of_mean_Tzero_err[i] = data[5]
        list_of_std_Tzero[i] = data[6]
        
        ######################################## prepare for hist of all beta_matrix_good        
        filename = 'beta_good_array_%s.dat' % (NP)
        beta_good_array_filepath = os.path.join(beta_folder, filename)
        data = np.loadtxt(beta_good_array_filepath)
        list_of_beta_good_array.append(data)
        key = 'col_%02d_np_%02d' % (list_of_col[i], list_of_nps[i])
        list_of_beta_good_dict[key] = data
        
        filename = 'beta_good_err_array_%s.dat' % (NP)
        beta_good_err_array_filepath = os.path.join(beta_folder, filename)
        data = np.loadtxt(beta_good_err_array_filepath)
        key = 'col_%02d_np_%02d' % (list_of_col[i], list_of_nps[i])
        list_of_beta_good_err_dict[key] = data
        
        ############################################ import beta
        Tzero_folder = os.path.join(path_from, NP, 'matrix', 'Tzero_matrix')       
        ######################################## prepare for hist of all Tzero_matrix_good        
        filename = 'Tzero_good_array_%s.dat' % (NP)
        Tzero_good_array_filepath = os.path.join(Tzero_folder, filename)
        data = np.loadtxt(Tzero_good_array_filepath)
        list_of_Tzero_good_array.append(data)
        key = 'col_%02d_np_%02d' % (list_of_col[i], list_of_nps[i])
        list_of_Tzero_good_dict[key] = data
        
        filename = 'Tzero_good_err_array_%s.dat' % (NP)
        Tzero_good_err_array_filepath = os.path.join(Tzero_folder, filename)
        data = np.loadtxt(Tzero_good_err_array_filepath)
        key = 'col_%02d_np_%02d' % (list_of_col[i], list_of_nps[i])
        list_of_Tzero_good_err_dict[key] = data
        
        ############################################ import spr data    
        spr_folder = os.path.join(path_from, NP, 'spr')
        filename = 'spr_fitted_parameters_%s.dat' % NP
        spr_filepath = os.path.join(spr_folder, filename)
        data = np.loadtxt(spr_filepath)
        list_of_londa_max[i] = data[0]
        list_of_width_pl[i] = data[1]
        list_of_r_sq_spr[i] = data[2]
        
        ############################################ import power monitor data    
        if monitor:
            monitor_folder = os.path.join(path_from, NP, 'power_monitor')
            filename = 'power_monitor_%s.dat' % NP
            monitor_filepath = os.path.join(monitor_folder, filename)
            data = np.loadtxt(monitor_filepath)
            list_of_avg_monitor[i] = data[0]
            list_of_std_monitor[i] = data[1]
    
    ########## PLOT
    ########## PLOT
    ########## PLOT
    
    plt.figure()
    nbins = 20
    range_tuple = [20,120]
#    print(list_of_beta_good_array)
    plt.hist(list_of_beta_good_array, bins=nbins, range=range_tuple, rwidth = 1, \
             align='mid', alpha = 1, edgecolor='k', normed = False,
             histtype = 'barstacked', stacked=True)
    plt.legend(loc=1, handlelength=0)
    ax = plt.gca()
#    ax.axvline(mu_slope, ymin = 0, ymax = 1, color='k', linestyle='--', linewidth=2)
#    ax.fill_between([borders_slope[0], borders_slope[1]], [0,0], [100,100],  facecolor='k', alpha=0.25)
    ax.get_yaxis().set_ticks([])
    ax.set_ylim([0, 200])
    ax.set_xlim(range_tuple)
    ax.set_xlabel(u'Photothermal coefficient (K µm$^{2}$/mW)')
    ax.set_ylabel('Entries')
    aux_folder = manage_save_directory(path_from,'stats')
    figure_name = os.path.join(aux_folder,'hist_of_all_betas_R2th_%s_%s_bines.png' % \
                               (str(R2th), str(totalbins)))
    plt.savefig(figure_name)
    plt.close()

    plt.figure()
    nbins = 20
    range_tuple = [250,350]
#    print(list_of_beta_good_array)
    plt.hist(list_of_Tzero_good_array, bins=nbins, range=range_tuple, rwidth = 1, \
             align='mid', alpha = 1, edgecolor='k', normed = False,
             histtype = 'barstacked', stacked=True)
    plt.legend(loc=1, handlelength=0)
    ax = plt.gca()
#    ax.axvline(mu_slope, ymin = 0, ymax = 1, color='k', linestyle='--', linewidth=2)
#    ax.fill_between([borders_slope[0], borders_slope[1]], [0,0], [100,100],  facecolor='k', alpha=0.25)
    ax.get_yaxis().set_ticks([])
    ax.set_ylim([0, 200])
    ax.set_xlim(range_tuple)
    ax.set_xlabel(u'T$_{0}$ (K)')
    ax.set_ylabel('Entries')
    aux_folder = manage_save_directory(path_from,'stats')
    figure_name = os.path.join(aux_folder,'hist_of_all_Tzero_R2th_%s_%s_bines.png' % \
                               (str(R2th), str(totalbins)))
    plt.savefig(figure_name)
    plt.close()

    #################### SAVE DATA ############################
    #################### SAVE DATA ############################
    #################### SAVE DATA ############################
    
    aux_folder = manage_save_directory(path_from,'stats')
    
    to_save = np.array([list_of_mean_beta,
                         list_of_mean_beta_err,
                         list_of_std_beta,
                         list_of_data_points,
                         list_of_mean_Tzero,
                         list_of_mean_Tzero_err,
                         list_of_std_Tzero,
                         list_of_londa_max, 
                         list_of_width_pl, 
                         list_of_r_sq_spr,
                         list_of_avg_monitor,
                         list_of_std_monitor,
                         list_of_col,
                         list_of_nps]).T
    header_text = 'BETA_mW_um2_K ERR_BETA_mW_um2_K STD_BETA_mW_um2_K DATA_POINTS TZERO_K ERR_TZERO_K STD_TZERO_K LAMBDA_MAX_nm WIDTH_nm SPR_R_SQ MONITOR_AVG MONITOR_STD COL NP'
    path_to_save = os.path.join(aux_folder,'all_NP_data.dat')
    np.savetxt(path_to_save, to_save, fmt='%.3e', header=header_text, comments='')
    
    path_to_save = os.path.join(aux_folder,'all_beta_good_data.pkl')
    f = open(path_to_save,'wb+')
    pickle.dump(list_of_beta_good_dict,f)
    f.close()
    
    path_to_save = os.path.join(aux_folder,'all_beta_good_err_data.pkl')
    f = open(path_to_save,'wb+')
    pickle.dump(list_of_beta_good_err_dict,f)
    f.close()

    path_to_save = os.path.join(aux_folder,'all_Tzero_good_data.pkl')
    f = open(path_to_save,'wb+')
    pickle.dump(list_of_Tzero_good_dict,f)
    f.close()
    
    path_to_save = os.path.join(aux_folder,'all_Tzero_good_err_data.pkl')
    f = open(path_to_save,'wb+')
    pickle.dump(list_of_Tzero_good_err_dict,f)
    f.close()
    
    return

def statistics(path_from, R2th, totalbins, radius, sigma_abs, find_beta, find_Tzero):
    
    stats_file = os.path.join(path_from,'stats','all_NP_data.dat')
    data = np.loadtxt(stats_file, skiprows=1)

    list_of_mean_beta = data[:,0]
    list_of_mean_beta_err = data[:,1]
    list_of_std_beta = data[:,2]
    list_of_data_points = data[:,3]
    list_of_mean_Tzero = data[:,4]
    list_of_mean_Tzero_err = data[:,5]
    list_of_std_Tzero = data[:,6]
    list_of_londa_max = data[:,7]
    list_of_width_pl = data[:,8]
    list_of_r_sq_spr = data[:,9]  
    list_of_avg_monitor = data[:,10]   
    list_of_std_monitor = data[:,11]
    list_of_col = data[:,12]   
    list_of_nps = data[:,13]
    
    index_non_zero = list_of_mean_beta!=0
    list_of_mean_beta = list_of_mean_beta[index_non_zero]
    list_of_mean_beta_err = list_of_mean_beta_err[index_non_zero]
    list_of_std_beta = list_of_std_beta[index_non_zero]
    list_of_data_points = list_of_data_points[index_non_zero]
    list_of_mean_Tzero = list_of_mean_Tzero[index_non_zero]
    list_of_mean_Tzero_err = list_of_mean_Tzero_err[index_non_zero]
    list_of_std_Tzero = list_of_std_Tzero[index_non_zero]
    list_of_londa_max = list_of_londa_max[index_non_zero]
    list_of_width_pl = list_of_width_pl[index_non_zero]
    list_of_r_sq_spr = list_of_r_sq_spr[index_non_zero]
    list_of_avg_monitor = list_of_avg_monitor[index_non_zero]
    list_of_std_monitor = list_of_std_monitor[index_non_zero]
    list_of_col = list_of_col[index_non_zero]
    list_of_nps = list_of_nps[index_non_zero]
   
    ##########################################################################
    ##########################################################################
    ##########################################################################
    
    ######################################## slope with errors vs NP
    plt.figure()
    plt.errorbar(list_of_nps, list_of_mean_beta, yerr=list_of_mean_beta_err, \
             elinewidth=1, fmt='s', alpha=0.5, label='All beta NPs')
    plt.legend()
    ax = plt.gca()
    ax.set_ylim([0, 100])
    ax.set_xlabel(r'Index')
    ax.set_ylabel(r'Photothermal coefficient (K µm$^{2}$/mW)')
    aux_folder = manage_save_directory(path_from,'stats')
    figure_name = 'beta_vs_index_R2th_%s_%s_bines.png' % \
                               (str(R2th), str(totalbins))
    figure_name = os.path.join(aux_folder, figure_name)
    plt.savefig(figure_name)
    plt.close()
    
    ######################################## slope vs NP
    plt.figure()
    plt.plot(list_of_col, list_of_mean_beta, 's', color='C1', label='All beta NPs')
    plt.legend()
    ax = plt.gca()
    ax.set_ylim([0, 100])
    ax.set_xlabel(r'NP')
    ax.set_ylabel(r'Photothermal coefficient (K µm$^{2}$/mW)')
    aux_folder = manage_save_directory(path_from,'stats')
    figure_name = 'beta_vs_NP_R2th_%s_%s_bines.png' % \
                               (str(R2th), str(totalbins))
    figure_name = os.path.join(aux_folder, figure_name)
    plt.savefig(figure_name)
    plt.close()
    
    ######################################## hist BETA 
    ######################################## hist BETA 
    ######################################## hist BETA 
        
    unique_NPs = np.unique(list_of_col)
    number_of_NPs = len(unique_NPs)
    range_tuple = [15,95]
    nbins = 16
    
    if find_beta and not find_Tzero:
        
        ####################################### hist BETA per NP
        #### ALL NPs IN ONE PLOT (histogram of its repetition)
        if number_of_NPs != 1:
            fig, axs = plt.subplots(number_of_NPs, 1, gridspec_kw = {'wspace':0, 'hspace':0})    
            for i in range(number_of_NPs):
                num_NP = int(unique_NPs[i])
                grep_NPs = np.where(list_of_col == num_NP)[0]
                grep_betas_NPs = list_of_mean_beta[grep_NPs]
                beta_mean_per_NP = np.nanmean(grep_betas_NPs)
                betas_std_per_NP = np.nanstd(grep_betas_NPs, ddof = 1)
                borders = [beta_mean_per_NP + betas_std_per_NP, beta_mean_per_NP - betas_std_per_NP]
                ax = axs[i]
                out_hist = ax.hist(grep_betas_NPs, bins=nbins, range=range_tuple, rwidth = 1, \
                     align='mid', alpha = 1, color = 'C0', edgecolor='k', normed = False, zorder=2, \
                     label='NP %d' % i)
                ax.axvline(beta_mean_per_NP, ymin = 0, ymax = 1, color='C3', linestyle='--', linewidth=2, zorder=3)
                ax.fill_between([borders[0], borders[1]], [0,0], [1,1],  facecolor='k', alpha=0.25, zorder=1)
                ax.set_xlim([41, 89])
            plt.legend()
            ax = plt.gca()
            ax.set_ylim([0, 10])
            ax.set_ylabel('Entries')
            ax.set_xlabel(u'Photothermal coefficient (K µm$^{2}$/mW)')
            aux_folder = manage_save_directory(path_from,'stats')
            figure_name = '0_beta_mean_hist_per_NP_R2th_%s_%s_bines.png' % \
                                       (str(R2th), str(totalbins))
            figure_name = os.path.join(aux_folder, figure_name)
            plt.savefig(figure_name, dpi = 300, bbox_layout = 'tight')
            plt.close()
        
            
        ####################################### hist BETA MEAN per SIZE
        #### ALL NPs IN ONE PLOT (histogram of its repetition)
        if number_of_NPs != 1:
            fig, axs = plt.subplots(number_of_NPs, 1, gridspec_kw = {'wspace':0, 'hspace':0})    
            for i in range(number_of_NPs):
                num_NP = int(unique_NPs[i])
                grep_NPs = np.where(list_of_col == num_NP)[0]
                grep_betas_NPs = list_of_mean_beta[grep_NPs]
                beta_mean_per_NP = np.nanmean(grep_betas_NPs)
                betas_std_per_NP = np.nanstd(grep_betas_NPs, ddof = 1)
                borders = [beta_mean_per_NP + betas_std_per_NP, beta_mean_per_NP - betas_std_per_NP]
                ax = axs[i]
                out_hist = ax.hist(grep_betas_NPs, bins=nbins, range=range_tuple, rwidth = 1, \
                     align='mid', alpha = 1, color = 'C0', edgecolor='k', normed = False, zorder=2, \
                     label='NP %d' % i)
                ax.axvline(beta_mean_per_NP, ymin = 0, ymax = 1, color='C3', linestyle='--', linewidth=2, zorder=3)
                ax.fill_between([borders[0], borders[1]], [0,0], [1,1],  facecolor='k', alpha=0.25, zorder=1)
                ax.set_xlim([41, 89])
            plt.legend()
            ax = plt.gca()
            ax.set_ylim([0, 10])
            ax.set_ylabel('Entries')
            ax.set_xlabel(u'Photothermal coefficient (K µm$^{2}$/mW)')
            aux_folder = manage_save_directory(path_from,'stats')
            figure_name = '0_beta_mean_hist_per_NP_R2th_%s_%s_bines.png' % \
                                       (str(R2th), str(totalbins))
            figure_name = os.path.join(aux_folder, figure_name)
            plt.savefig(figure_name, dpi = 300, bbox_layout = 'tight')
            plt.close()
        
        
        ####################################### hist BETA per NP
        #### one NP per plot (histogram of its repetition)
        
        print('\n-------- Stats using matrix betas of multiple scans of one NP (to find distributions parameters) --------')
        
        dict_file = os.path.join(path_from,'stats','all_beta_good_data.pkl')
        with open(dict_file, 'rb') as f:
            list_of_beta_good_dict = pickle.load(f)
        dict_file = os.path.join(path_from,'stats','all_beta_good_err_data.pkl')
        with open(dict_file, 'rb') as f:
            list_of_beta_good_err_dict = pickle.load(f)
        for i in range(number_of_NPs):
            list_of_beta_good_per_NP = np.array([])
            list_of_beta_good_err_per_NP = np.array([])
            num_NP = int(unique_NPs[i])
            for key in list_of_beta_good_dict.keys():
                col = int(key.split('_')[1])
    #            NP = int(key.split('_')[3])
                if num_NP == col:
    #                print(num_NP, col, NP)
                    list_of_beta_good_per_NP = np.append(list_of_beta_good_dict[key], list_of_beta_good_per_NP)
                    list_of_beta_good_err_per_NP = np.append(list_of_beta_good_err_dict[key], list_of_beta_good_err_per_NP)
    #        print(list_of_beta_good_per_NP)
            
            # stats
            data_points = len(list_of_beta_good_per_NP)
            beta_mean_per_NP = np.nanmean(list_of_beta_good_per_NP)
            betas_std_per_NP = np.nanstd(list_of_beta_good_per_NP, ddof = 1)
    #        beta_mean_err_per_NP = np.sqrt(np.nansum(list_of_beta_good_err_per_NP**2))/data_points
            beta_mean_err_per_NP = betas_std_per_NP/np.sqrt(data_points)
            borders = [beta_mean_per_NP - betas_std_per_NP, beta_mean_per_NP + betas_std_per_NP]
            
            print('\nNP %02d' % num_NP)
            print('Data points (all betas of same NP scanned multiple times): %d' % data_points)
            print('Mean beta: %.3f' % beta_mean_per_NP, 'Std dev beta: %.3f' % betas_std_per_NP)
            print('Error in mean beta (same sigma for all): %.3f' % beta_mean_err_per_NP)        
            
            plt.figure()
            ax = plt.gca()
            out_hist = plt.hist(list_of_beta_good_per_NP, bins=nbins, range=range_tuple, rwidth = 1, \
                 align='mid', alpha = 1, color = 'C0', edgecolor='k', normed = False, zorder=2, \
                 label='NP %d' % i)
            ax.axvline(beta_mean_per_NP, ymin = 0, ymax = 1, color='C3', linestyle='--', linewidth=2, zorder=3)
            ax.fill_between([borders[0], borders[1]], [0,0], [100,100],  facecolor='k', alpha=0.25, zorder=1)
            plt.legend()
            ax.set_ylim([0, 1.1*np.max(out_hist[0])])
            ax.set_xlim([25, 105])
            ax.set_ylabel('Entries')
            ax.set_xlabel(r'$\beta$ (K µm$^{2}$/mW)')
            aux_folder = manage_save_directory(path_from,'stats')
            figure_name = '0_beta_hist_NP_%02d_R2th_%s_%s_bines.png' % \
                                       (num_NP, str(R2th), str(totalbins))
            figure_name = os.path.join(aux_folder, figure_name)
            plt.savefig(figure_name, dpi = 300, bbox_layout = 'tight')
            plt.close()
        
        
        print('\n-------- Stats using mean betas for one NP (repetition, to find beta for each NP) --------')
        
    #    range_tuple = [15,85]
    #    nbins = 15
        
        list_of_beta_bb = []
        list_of_beta_barra_std = []
        list_of_beta_bb_err = []
        for i in range(number_of_NPs):
            num_NP = int(unique_NPs[i])
            
            # keep NP data (usually is all daata from one Col)
            grep_NPs = np.where(list_of_col == num_NP)[0]
            grep_betas_NPs = list_of_mean_beta[grep_NPs]
            grep_betas_err_NPs = list_of_mean_beta_err[grep_NPs]
            
            # stats
            data_points = len(grep_betas_NPs)
            beta_barra_mean_per_NP = np.nanmean(grep_betas_NPs)
            beta_barra_std_per_NP = np.nanstd(grep_betas_NPs, ddof = 1)
            beta_barra_mean_err_per_NP = np.sqrt(np.nansum(grep_betas_err_NPs**2))/data_points
            borders = [beta_barra_mean_per_NP - beta_barra_std_per_NP, beta_barra_mean_per_NP + beta_barra_std_per_NP]
            
            print('\nNP %02d' % num_NP)
            print('Data points (number of scans): %d' % data_points)
            print('Mean beta barra (beta bb): %.3f' % beta_barra_mean_per_NP, 'Std dev beta barra: %.3f' % beta_barra_std_per_NP)
            print('Error in mean beta barra (same sigma): %.3f' % beta_barra_mean_err_per_NP)        
            
            list_of_beta_bb.append(beta_barra_mean_per_NP)
            list_of_beta_barra_std.append(beta_barra_std_per_NP)
            list_of_beta_bb_err.append(beta_barra_mean_err_per_NP)
            
            plt.figure()
            ax = plt.gca()
            out_hist = plt.hist(grep_betas_NPs, bins=nbins, range=range_tuple, rwidth = 1, \
                 align='mid', alpha = 1, color = 'C0', edgecolor='k', normed = False, zorder=2, \
                 label='NP %d' % i)
            ax.axvline(beta_barra_mean_per_NP, ymin = 0, ymax = 1, color='C3', linestyle='--', linewidth=2, zorder=3)
            ax.fill_between([borders[0], borders[1]], [0,0], [100,100],  facecolor='k', alpha=0.25, zorder=1)
            plt.legend()
            ax.set_ylim([0, 1.1*np.max(out_hist[0])])
            ax.set_xlim([39, 91])
            ax.set_ylabel('Entries')
            ax.set_xlabel(r'$\bar{\beta}$ (K µm$^{2}$/mW)')
            aux_folder = manage_save_directory(path_from,'stats')
            figure_name = '0_beta_mean_hist_NP_%02d_R2th_%s_%s_bines.png' % \
                                       (num_NP, str(R2th), str(totalbins))
            figure_name = os.path.join(aux_folder, figure_name)
            plt.savefig(figure_name, dpi = 300, bbox_layout = 'tight')
            plt.close()
            
    #    # FOR INSET
    ##    range_tuple = [15,85]
    ##    nbins = 15
    #    for i in range(number_of_NPs):
    #        num_NP = int(unique_NPs[i])
    #        grep_NPs = np.where(list_of_col == num_NP)[0]
    #        grep_betas_NPs = list_of_mean_beta[grep_NPs]
    #        beta_mean_per_NP = np.nanmean(grep_betas_NPs)
    #        betas_std_per_NP = np.nanstd(grep_betas_NPs, ddof = 1)
    #        borders = [beta_mean_per_NP - betas_std_per_NP, beta_mean_per_NP + betas_std_per_NP]
    #        plt.figure()
    #        ax = plt.gca()
    #        out_hist = plt.hist(grep_betas_NPs, bins=nbins, range=range_tuple, rwidth = 1, \
    #             align='mid', alpha = 1, color = 'C0', edgecolor='k', normed = False, zorder=2, \
    #             label='NP %d' % i)
    #        ax.axvline(beta_mean_per_NP, ymin = 0, ymax = 1, color='C3', linestyle='--', linewidth=2, zorder=3)
    #        ax.fill_between([borders[0], borders[1]], [0,0], [100,100],  facecolor='k', alpha=0.25, zorder=1)
    #        plt.legend()
    #        ax.set_ylim([0, 1.1*np.max(out_hist[0])])
    #        ax.set_xlim([49, 81])
    #        ax.tick_params(axis='both', which='major', labelsize=25)
    #        ax.set_ylabel('Entries', fontsize=24)
    #        ax.set_xlabel(r'$\bar{beta}$ (K µm$^{2}$/mW)', fontsize=24)
    #        aux_folder = manage_save_directory(path_from,'stats')
    #        fig_name = '0_beta_mean_hist_NP_%02d_R2th_%s_%s_bines_%s_Tzero_PAPER_INSET' % \
    #                                   (num_NP, str(R2th), str(totalbins), str(Tzero))
    #        figurefull_name = '%s.png' % fig_name
    #        figure_name = os.path.join(aux_folder, figurefull_name)
    #        plt.savefig(figure_name, dpi = 300, bbox_layout = 'tight')
    #        figurefull_name = '%s.pdf' % fig_name
    #        figure_name = os.path.join(aux_folder, figurefull_name)
    #        plt.savefig(figure_name, dpi = 300, bbox_layout = 'tight', format = 'pdf')
    #        plt.close()
          
          
        ######################################## slope vs max lambda
        plt.figure()
        plt.plot(list_of_londa_max, list_of_mean_beta, 'o', alpha=0.5, label='All beta NPs')
        plt.legend()
        ax = plt.gca()
        ax.set_ylim([0, 100])
        ax.set_xlim([530, 600])
        ax.set_xlabel(r'Max PL wavelength (nm)')
        ax.set_ylabel(r'Photothermal coefficient (K µm$^{2}$/mW)')
        aux_folder = manage_save_directory(path_from,'stats')
        figure_name = 'beta_vs_lambda_R2th_%s_%s_bines.png' % \
                                   (str(R2th), str(totalbins))
        figure_name = os.path.join(aux_folder, figure_name)
        plt.savefig(figure_name)
        plt.close()
                
        ######################################## power monitor vs slope, fake correlation
        plt.figure()
        plt.errorbar(list_of_avg_monitor, list_of_mean_beta, xerr = list_of_std_monitor, \
                 yerr = list_of_mean_beta_err, elinewidth=1, fmt='o', alpha=0.5, label='All NPs')
        plt.legend()
        ax = plt.gca()
        ax.set_ylim([0, 110])
        ax.set_ylabel(r'Photothermal coefficient (K µm$^{2}$/mW)')
        ax.set_xlabel(r'Power monitor (a.u.)')
        aux_folder = manage_save_directory(path_from,'stats')
        figure_name = os.path.join(aux_folder,'power_monitor_vs_slope_R2th_%s_%s_bines.png' % \
                                   (str(R2th), str(totalbins)))
        plt.savefig(figure_name)
        plt.close()

        print('\n-------- Stats using ALL beta barra bracket (one per NP, to find beta for a given system NP+medium) --------')
        
        number_of_beta_bb = len(list_of_beta_bb)
        print('Number of NPs measured: %d' % number_of_beta_bb)
        
        mu_beta_bb = np.mean(list_of_beta_bb)
        std_beta_bb = np.std(list_of_beta_bb, ddof=1)
        print('\nMean beta bb: %.3f' % mu_beta_bb, 'Std dev beta bb: %.3f' % std_beta_bb)

#        err_mu_beta_bb = np.sqrt(np.sum(beta_barra_mean_err_per_NP**2))/number_of_beta_bb
#        err_mu_beta_bb = np.sqrt(number_of_beta_bb)
#        beta_barra_std_per_NP
#        print('Error in mean beta bb: %.3f' % err_mu_slope)
        
        k_eff = sigma_abs/(4*np.pi*radius*mu_beta_bb)
        print('\nEffective thermal conductivity: %.3f' % k_eff, 'W/Km')

        

    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
   
    elif not find_beta and find_Tzero:
    
        ######################################## Tzero with errors vs NP
        plt.figure()
        plt.errorbar(list_of_nps, list_of_mean_beta, yerr=list_of_mean_beta_err, \
                 elinewidth=1, fmt='s', alpha=0.5, label='All Tzero NPs')
        plt.legend()
        ax = plt.gca()
        ax.set_ylim([0, 100])
        ax.set_xlabel(r'Index')
        ax.set_ylabel(r'T$_{0}$ (K)')
        aux_folder = manage_save_directory(path_from,'stats')
        figure_name = 'Tzero_vs_index_R2th_%s_%s_bines.png' % \
                                   (str(R2th), str(totalbins))
        figure_name = os.path.join(aux_folder, figure_name)
        plt.savefig(figure_name)
        plt.close()
        
        ######################################## Tzero vs NP
        plt.figure()
        plt.plot(list_of_col, list_of_mean_Tzero, 's', color='C1', label='All Tzero NPs')
        plt.legend()
        ax = plt.gca()
        ax.set_ylim([0, 100])
        ax.set_xlabel(r'NP')
        ax.set_ylabel(r'T$_{0}$ (K)')
        aux_folder = manage_save_directory(path_from,'stats')
        figure_name = 'Tzero_vs_NP_R2th_%s_%s_bines.png' % \
                                   (str(R2th), str(totalbins))
        figure_name = os.path.join(aux_folder, figure_name)
        plt.savefig(figure_name)
        plt.close()
        
        ######################################## hist Tzero per NP
        #### ALL NPs IN ONE (histogram of its repetition)
        
        unique_NPs = np.unique(list_of_col)
        number_of_NPs = len(unique_NPs)
        range_tuple = [280,380]
        nbins = 16
        if number_of_NPs != 1:
            fig, axs = plt.subplots(number_of_NPs, 1, gridspec_kw = {'wspace':0, 'hspace':0})    
            for i in range(number_of_NPs):
                num_NP = int(unique_NPs[i])
                grep_NPs = np.where(list_of_col == num_NP)[0]
                grep_Tzero_NPs = list_of_mean_Tzero[grep_NPs]
                Tzero_mean_per_NP = np.mean(grep_Tzero_NPs)
                Tzero_std_per_NP = np.std(grep_Tzero_NPs, ddof = 1)
                borders = [Tzero_mean_per_NP + Tzero_std_per_NP, Tzero_mean_per_NP - Tzero_std_per_NP]
                ax = axs[i]
                out_hist = ax.hist(grep_Tzero_NPs, bins=nbins, range=range_tuple, rwidth = 1, \
                     align='mid', alpha = 1, color = 'C0', edgecolor='k', normed = False, zorder=2, \
                     label='NP %d' % i)
                ax.axvline(Tzero_mean_per_NP, ymin = 0, ymax = 1, color='C3', linestyle='--', linewidth=2, zorder=3)
                ax.fill_between([borders[0], borders[1]], [0,0], [1,1],  facecolor='k', alpha=0.25, zorder=1)
                ax.set_xlim([280, 380])
            plt.legend()
            ax = plt.gca()
            ax.set_ylim([0, 1.1*np.max(out_hist[0])])
            ax.set_ylabel('Entries')
            ax.set_xlabel(u'T$_{0}$ (K)')
            aux_folder = manage_save_directory(path_from,'stats')
            figure_name = 'Tzero_hist_per_NP_R2th_%s_%s_bines.png' % \
                                       (str(R2th), str(totalbins))
            figure_name = os.path.join(aux_folder, figure_name)
            plt.savefig(figure_name, dpi = 300, bbox_layout = 'tight')
            plt.close()
        
        ####################################### hist BETA per NP
        #### one NP per plot (histogram of its repetition)
        
        print('\n-------- Stats using matrix Tzero (not mean) for one NP --------')
        
        dict_file = os.path.join(path_from,'stats','all_Tzero_good_data.pkl')
        with open(dict_file, 'rb') as f:
            list_of_Tzero_good_dict = pickle.load(f)
        dict_file = os.path.join(path_from,'stats','all_Tzero_good_err_data.pkl')
        with open(dict_file, 'rb') as f:
            list_of_Tzero_good_err_dict = pickle.load(f)
        for i in range(number_of_NPs):
            list_of_Tzero_good_per_NP = np.array([])
            list_of_Tzero_good_err_per_NP = np.array([])
            num_NP = int(unique_NPs[i])
            for key in list_of_Tzero_good_dict.keys():
                col = int(key.split('_')[1])
    #            NP = int(key.split('_')[3])
                if num_NP == col:
    #                print(num_NP, col, NP)
                    list_of_Tzero_good_per_NP = np.append(list_of_Tzero_good_dict[key], list_of_Tzero_good_per_NP)
                    list_of_Tzero_good_err_per_NP = np.append(list_of_Tzero_good_err_dict[key], list_of_Tzero_good_err_per_NP)
    #        print(list_of_Tzero_good_per_NP)
            
            # stats
            data_points = len(list_of_Tzero_good_per_NP)
            Tzero_mean_per_NP = np.mean(list_of_Tzero_good_per_NP)
            Tzero_std_per_NP = np.std(list_of_Tzero_good_per_NP, ddof = 1)
    #        Tzero_mean_err_per_NP = np.sqrt(np.sum(list_of_Tzero_good_err_per_NP**2))/data_points
            Tzero_mean_err_per_NP = Tzero_std_per_NP/np.sqrt(data_points)
            borders = [Tzero_mean_per_NP - Tzero_std_per_NP, Tzero_mean_per_NP + Tzero_std_per_NP]
            
            print('\nNP %02d' % num_NP)
            print('Data points: %d' % data_points)
            print('Mean Tzero: %.3f' % Tzero_mean_per_NP, 'Std dev Tzero: %.3f' % Tzero_std_per_NP)
            print('Error in mean Tzero (same sigma for all): %.3f' % Tzero_mean_err_per_NP)        
            
            plt.figure()
            ax = plt.gca()
            out_hist = plt.hist(list_of_Tzero_good_per_NP, bins=nbins, range=range_tuple, rwidth = 1, \
                 align='mid', alpha = 1, color = 'C0', edgecolor='k', normed = False, zorder=2, \
                 label='NP %d' % i)
            ax.axvline(Tzero_mean_per_NP, ymin = 0, ymax = 1, color='C3', linestyle='--', linewidth=2, zorder=3)
            ax.fill_between([borders[0], borders[1]], [0,0], [100,100],  facecolor='k', alpha=0.25, zorder=1)
            plt.legend()
            ax.set_ylim([0, 1.1*np.max(out_hist[0])])
            ax.set_xlim([280, 380])
            ax.set_ylabel('Entries')
            ax.set_xlabel(u'T$_{0}$ (K)')
            aux_folder = manage_save_directory(path_from,'stats')
            figure_name = 'Tzero_hist_NP_%02d_R2th_%s_%s_bines.png' % \
                                       (num_NP, str(R2th), str(totalbins))
            figure_name = os.path.join(aux_folder, figure_name)
            plt.savefig(figure_name, dpi = 300, bbox_layout = 'tight')
            plt.close()
        
        
        print('\n-------- Stats using mean Tzero for one NP --------')
        
    #    range_tuple = [15,85]
    #    nbins = 15
        list_of_Tzero_bb = []
        list_of_Tzero_barra_std = []
        list_of_Tzero_bb_err = []
        for i in range(number_of_NPs):
            num_NP = int(unique_NPs[i])
            
            # keep NP data (usually is all daata from one Col)
            grep_NPs = np.where(list_of_col == num_NP)[0]
            grep_Tzero_NPs = list_of_mean_Tzero[grep_NPs]
            grep_Tzero_err_NPs = list_of_mean_Tzero_err[grep_NPs]
            
            # stats
            data_points = len(grep_Tzero_NPs)
            Tzero_barra_mean_per_NP = np.mean(grep_Tzero_NPs)
            Tzero_barra_std_per_NP = np.std(grep_Tzero_NPs, ddof = 1)
            Tzero_barra_mean_err_per_NP = np.sqrt(np.sum(grep_Tzero_err_NPs**2))/data_points
            borders = [Tzero_barra_mean_per_NP - Tzero_barra_std_per_NP, Tzero_barra_mean_per_NP + Tzero_barra_std_per_NP]
            
            print('\nNP %02d' % num_NP)
            print('Data points: %d' % data_points)
            print('Mean Tzero barra (Tzero bb): %.3f' % Tzero_barra_mean_per_NP, 'Std dev Tzero barra: %.3f' % Tzero_barra_std_per_NP)
            print('Error in mean Tzero barra (not same sigma): %.3f' % Tzero_barra_mean_err_per_NP)        
            
            list_of_Tzero_bb.append(Tzero_barra_mean_per_NP)
            list_of_Tzero_barra_std.append(Tzero_barra_std_per_NP)
            list_of_Tzero_bb_err.append(Tzero_barra_mean_err_per_NP)
            
            plt.figure()
            ax = plt.gca()
            out_hist = plt.hist(grep_Tzero_NPs, bins=nbins, range=range_tuple, rwidth = 1, \
                 align='mid', alpha = 1, color = 'C0', edgecolor='k', normed = False, zorder=2, \
                 label='NP %d' % i)
            ax.axvline(Tzero_barra_mean_per_NP, ymin = 0, ymax = 1, color='C3', linestyle='--', linewidth=2, zorder=3)
            ax.fill_between([borders[0], borders[1]], [0,0], [100,100],  facecolor='k', alpha=0.25, zorder=1)
            plt.legend()
            ax.set_ylim([0, 1.1*np.max(out_hist[0])])
            ax.set_xlim([280, 380])
            ax.set_ylabel('Entries')
            ax.set_xlabel(u'T$_{0}$ (K)')
            aux_folder = manage_save_directory(path_from,'stats')
            figure_name = 'Tzero_mean_hist_NP_%02d_R2th_%s_%s_bines.png' % \
                                       (num_NP, str(R2th), str(totalbins))
            figure_name = os.path.join(aux_folder, figure_name)
            plt.savefig(figure_name, dpi = 300, bbox_layout = 'tight')
            plt.close()
          
        ######################################## Tzero vs max lambda
        plt.figure()
        plt.plot(list_of_londa_max, list_of_mean_Tzero, 'o', alpha=0.5, label='All Tzero NPs')
        plt.legend()
        ax = plt.gca()
        ax.set_ylim([250, 400])
        ax.set_xlim([530, 600])
        ax.set_xlabel(r'Max PL wavelength (nm)')
        ax.set_ylabel(r'T$_{0}$ (K)')
        aux_folder = manage_save_directory(path_from,'stats')
        figure_name = 'Tzero_vs_lambda_R2th_%s_%s_bines.png' % \
                                   (str(R2th), str(totalbins))
        figure_name = os.path.join(aux_folder, figure_name)
        plt.savefig(figure_name)
        plt.close()
        
        ######################################## power monitor vs Tzero, fake correlation
        plt.figure()
        plt.errorbar(list_of_avg_monitor, list_of_mean_Tzero, xerr = list_of_std_monitor, \
                 yerr = list_of_mean_Tzero_err, elinewidth=1, fmt='o', alpha=0.5, label='All NPs')
        plt.legend()
        ax = plt.gca()
        ax.set_ylim([0, 110])
        ax.set_ylabel(r'T$_{0}$ (K)')
        ax.set_xlabel(r'Power monitor (a.u.)')
        aux_folder = manage_save_directory(path_from,'stats')
        figure_name = os.path.join(aux_folder,'power_monitor_vs_Tzero_R2th_%s_%s_bines.png' % \
                                   (str(R2th), str(totalbins)))
        plt.savefig(figure_name)
        plt.close()
        
        print('\n-------- Stats using ALL Tzero barra bracket (one per NP, to find Tzero for a given system NP+medium) --------')
        
        number_of_Tzero_bb = len(list_of_Tzero_bb)
        print('Number of NPs measured: %d' % number_of_Tzero_bb)
        
        mu_Tzero_bb = np.mean(list_of_Tzero_bb)
        std_Tzero_bb = np.std(list_of_Tzero_bb, ddof=1)
        print('\nMean Tzero bb: %.3f' % mu_Tzero_bb, 'Std dev Tzero bb: %.3f' % std_Tzero_bb)

#        err_mu_beta_bb = np.sqrt(np.sum(beta_barra_mean_err_per_NP**2))/number_of_beta_bb
#        err_mu_beta_bb = np.sqrt(number_of_beta_bb)
#        beta_barra_std_per_NP
#        print('Error in mean beta bb: %.3f' % err_mu_slope)        

    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################    
    
    ######################################## lambda max vs NP
    plt.figure()
    plt.plot(list_of_nps, list_of_londa_max, 'o')
    plt.legend()
    ax = plt.gca()
    ax.set_ylim([530, 600])
    ax.set_xlabel(r'Index')
    ax.set_ylabel(r'Max PL wavelength (nm)')
    ax.axvline(9.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(19.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(29.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(39.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(49.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(59.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(69.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(79.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(89.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    aux_folder = manage_save_directory(path_from,'stats')
    figure_name = os.path.join(aux_folder,'max_lambda_vs_index_R2th_%s_%s_bines.png' % \
                               (str(R2th), str(totalbins)))
    plt.savefig(figure_name)
    plt.close()

    ######################################## power monitor vs NP
    plt.figure()
    plt.errorbar(list_of_nps, list_of_avg_monitor, yerr = list_of_std_monitor, \
             elinewidth=1, fmt='o', alpha=0.5, label='All NPs')
    plt.legend()
    ax = plt.gca()
    ax.set_ylim([0.1, 0.2])
    ax.set_xlabel(r'Index')
    ax.set_ylabel(r'Power monitor (a.u.)')
    ax.axvline(9.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(19.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(29.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(39.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(49.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(59.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(69.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(79.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    ax.axvline(89.5, ymin=0, ymax=1, color='k', alpha=0.5, linestyle='--')
    aux_folder = manage_save_directory(path_from,'stats')
    figure_name = os.path.join(aux_folder,'power_monitor_vs_index_R2th_%s_%s_bines.png' % \
                               (str(R2th), str(totalbins)))
    plt.savefig(figure_name)
    plt.close()
    
    ########################################
    ########################################
    ########################################
    
    return