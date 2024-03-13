#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 09:53:36 2024

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
import sys
import copy, matplotlib
sys.path.insert(0, '../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info

prop_name = 'magnitude'
prop_name = 'R_sersic'
# prop_name = 'n_sersic'
psf_sp = 1 #!!! for the new run
all_library = True #!!! for the new run that use all the PSF library, or just only in FOV.
# add_cond = '_fixn1'  #!!!
add_cond = ''  #!!!

best_chisq_list_2022 = []
best_chisq_list = []
plt.figure(figsize=(11.5,8))
for idx in range(10):
    info = target_info[str(idx)]
    target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
    run_folder_2022 = '../../2022_JWST_QSOz6/model_z6_data_id{0}/stage3_all/'.format(idx) #!!!
    z_str = str(z)
    filt = ['F356W'][0]
    #%%Get first run result
    top_psf_id = 0
    first_fit_run_list = []
    if filt == 'F150W' :
        cmap = 'inferno'
    else:
        cmap = 'gist_heat'
    my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap)) # copy the default cmap
    my_cmap.set_bad('black')
    PSF_lib_files = glob.glob(run_folder_2022+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
    if idx ==1:
        fit_files = glob.glob(run_folder_2022+'*fit_material/fit_run*_fixn1_*idx{0}_{1}_*.pkl'.format(idx, filt))#+\
    else:
        fit_files = glob.glob(run_folder_2022+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
    fit_files.sort()
    
    for i in range(len(fit_files)):
        first_fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([first_fit_run_list[i].reduced_Chisq for i in range(len(first_fit_run_list))])
    sort_Chisq = chisqs.argsort()  
    fit_run = first_fit_run_list[sort_Chisq[top_psf_id]]
    count_n = 5
    Chisq_best = chisqs[sort_Chisq[top_psf_id]]
    Chisq_last= chisqs[sort_Chisq[count_n-1]]
    inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
    weight = np.zeros(len(chisqs))
    for i in sort_Chisq[:count_n]:
        weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
    all_values = [first_fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(first_fit_run_list))]
    weighted_value_2022 = np.sum(np.array(all_values)*weight) / np.sum(weight)
    best_chisq_list_2022.append(fit_run.reduced_Chisq)
    # rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value_2022)**2*weight) / np.sum(weight))
    if idx == 0:
        fit_run.plot_final_qso_fit(target_ID = target_id+'$-$'+filt, cmap = my_cmap)
        print(fit_run.final_result_galaxy[0], fit_run.reduced_Chisq)
    
    #%%Get second run result
    run_folder = '../material/fit_result/'
    if all_library == True:
        load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx*_psfsf{2}{3}.pkl'.format(filt, idx, psf_sp, add_cond))
    else:
        psf_idx = idx
        load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx{2}*_psfsf{3}{4}.pkl'.format(filt, idx, psf_idx, psf_sp,add_cond))
    load_files.sort()
    chisqs_idx = []
    fit_run_list = []
    for file in load_files:
        fit_run_list.append(pickle.load(open(file,'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort()  
    fit_run = fit_run_list[sort_Chisq[0]]
    count_n = 5
    all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
    weight = np.zeros(len(chisqs))
    for i in sort_Chisq[:count_n]:
        weight[i] = 1
    weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
    rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
    best_chisq_list.append(fit_run.reduced_Chisq)
    # if idx == 0:
    #     fit_run.plot_final_qso_fit(target_ID = target_id+'$-$'+filt, cmap = my_cmap)
    #     print(fit_run.final_result_galaxy[0], fit_run.reduced_Chisq)
    
    if prop_name == 'magnitude':  
        plt.scatter(idx, weighted_value-weighted_value_2022, c = fit_run.final_result_galaxy[0]['magnitude'], vmin=20, vmax=30)
        plt.errorbar(idx, weighted_value-weighted_value_2022, yerr=rms_value, ecolor = 'black')
        plt.ylim(-4,4)
    if prop_name == 'R_sersic':
        plt.scatter(idx, np.log10(weighted_value)-np.log10(weighted_value_2022), c = fit_run.final_result_galaxy[0]['magnitude'], vmin=20, vmax=30)
        plt.errorbar(idx, np.log10(weighted_value)-np.log10(weighted_value_2022), yerr=np.log10(weighted_value+rms_value) - np.log10(weighted_value), ecolor = 'black')
        plt.ylim(-0.5,0.5)
#%%
if prop_name == 'R_sersic':
    prop_name = '(in log)' + prop_name 
    plt.axhline(y=-0.2, linestyle ='--', linewidth = 2)
    plt.axhline(y=0.2, linestyle ='--', linewidth = 2)
if prop_name == 'magnitude':  
    plt.axhline(y=-0.4, linestyle ='--', linewidth = 2)
    plt.axhline(y=0.4, linestyle ='--', linewidth = 2)
plt.xlabel("Target ID",fontsize=27)
plt.ylabel("{0} (new-previous)".format(prop_name),fontsize=27)
cbar = plt.colorbar()
#cbar.set_clim(-2.0, 2.0)
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_ylabel('Host mag', rotation=270, fontsize = 25, labelpad=25)
plt.tick_params(labelsize=20)
plt.show()

if all_library == False:
    best_chisq_list_onlyFOV = best_chisq_list #!!! Run with all_library = False
#%%
plt.figure(figsize=(11.5,8))
for i in range(len(best_chisq_list)):
    plt.scatter(i, (best_chisq_list_onlyFOV[i] - best_chisq_list[i]) / best_chisq_list_onlyFOV[i] * 100 , c = fit_run.final_result_galaxy[0]['magnitude'], vmin=20, vmax=30)
plt.xlabel("Target ID",fontsize=27)
plt.ylabel("Improvement of chisq (%)".format(prop_name),fontsize=27)
plt.tick_params(labelsize=20)
plt.show()