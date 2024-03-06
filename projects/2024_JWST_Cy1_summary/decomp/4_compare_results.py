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

idx = 0
plt.figure(figsize=(11.5,8))
for idx in range(10):
    # prop_name = 'magnitude'
    # prop_name = 'R_sersic'
    prop_name = 'n_sersic'
    info = target_info[str(idx)]
    target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
    run_folder_2022 = '../../2022_JWST_QSOz6/model_z6_data_id{0}/stage3_all/'.format(idx) #!!!
    z_str = str(z)
    filt = ['F356W'][0]
    #%%Get first run result
    top_psf_id = 0
    fit_run_list = []
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
        fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort()  
    # print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
    fit_run = fit_run_list[sort_Chisq[top_psf_id]]
    # fit_run.plot_final_qso_fit(target_ID = target_id+'$-$'+filt, save_plot = False, cmap = my_cmap)
    count_n = 5
    Chisq_best = chisqs[sort_Chisq[top_psf_id]]
    Chisq_last= chisqs[sort_Chisq[count_n-1]]
    inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
    weight = np.zeros(len(chisqs))
    for i in sort_Chisq[:count_n]:
        weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
    # if idx == 2:
    #     fit_run.plot_final_qso_fit(target_ID = target_id+'$-$'+filt, cmap = my_cmap)
    all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
    weighted_value_2022 = np.sum(np.array(all_values)*weight) / np.sum(weight)
    rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value_2022)**2*weight) / np.sum(weight))
    host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
    AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
    ratio = host_flux/(host_flux+AGN_flux)
    # print(weighted_value_2022)
    
    #%%Get second run result
    run_folder = '../material/fit_result/'
    psf_sp = 2 #!!!
    # all_library = False
    all_library = True
    psf_idx = idx
    check_files = []
    add_cond = ''
    if idx == 1:
        add_cond = '_fixn1' 
    if all_library == True:
        load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx*_psfsf{2}{3}.pkl'.format(filt, idx, psf_sp, add_cond))
    else:
        load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx{2}*_psfsf{3}{4}.pkl'.format(filt, idx, psf_idx, psf_sp,add_cond))
    load_files.sort()
    chisqs_idx = []
    for file in load_files:
        fit_run_list.append(pickle.load(open(file,'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort()  
    fit_run = fit_run_list[sort_Chisq[0]]
    count_n = 3
    all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
    weight = np.zeros(len(chisqs))
    for i in sort_Chisq[:count_n]:
        weight[i] = 1
    weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
    rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
    plt.scatter(idx, weighted_value-weighted_value_2022, c = fit_run.final_result_galaxy[0]['magnitude'], vmin=20, vmax=30)
    plt.errorbar(idx, weighted_value-weighted_value_2022, yerr=rms_value, ecolor = 'black')

plt.ylim(-4,4)
plt.xlabel("Target ID",fontsize=27)
plt.ylabel("{0} (new-previous)".format(prop_name),fontsize=27)
cbar = plt.colorbar()
#cbar.set_clim(-2.0, 2.0)
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_ylabel('Host mag', rotation=270, fontsize = 25, labelpad=25)
plt.tick_params(labelsize=20)
plt.show()
