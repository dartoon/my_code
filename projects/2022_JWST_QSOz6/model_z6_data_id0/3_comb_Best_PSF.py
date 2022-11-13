#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 15:03:04 2022

@author: Dartoon
"""


import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle
import copy


run_folder = 'stage3_all_largecut/' #!!!
filt = 'F356W'
files = glob.glob(run_folder+'fit_material/data_process_idx0_*{0}*_*FOVpsf*.pkl'.format(filt))
files.sort()
collect_info = []
for i in range(len(files)):
    _file = files[i]
    idx_info = _file.split('idx')[1].split('_')[0]
    filt_info = _file.split('W_')[0].split('_')[-1] + 'W'
    this_info = [idx_info, filt_info]
    if this_info not in collect_info:
        collect_info.append(this_info)

# - [ ] F115W_psf6 is QSO (idx 2). [136, 137, 139]
# - [ ] F150W_psf7 is QSO (idx 2).
# - [ ] F277W_psf2 is QSO (idx 2).

#%%
if_printshow = False
for count in range(len(collect_info)):
    item = collect_info[count]
    fit_run_list = []
    idx, filt= item
    fit_files = glob.glob(run_folder+'fit_material/fit_run_idx{0}_{1}_*FOVpsf*.pkl'.format(idx, filt))
    fit_files.sort()
    
    for i in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    idx_counts = chisqs.argsort()  
    if len(idx_counts)<8:
        print(idx, filt, len(idx_counts))
    print("work on", count, 'idx', idx, filt, "Total PSF NO.", len(idx_counts))
    for ct in [2, 3, 5]:
        _data_process_list = [fit_run_list[i].fitting_specify_class.data_process_class for i in idx_counts[:ct]] 
        PSF_list_for_comb = [_data_process_list[i].PSF_list[0] for i in range(len(_data_process_list))]
        _data_process = copy.deepcopy(_data_process_list[0])
        _data_process.PSF_list  = copy.deepcopy(PSF_list_for_comb)
        _data_process.stack_PSF(if_plot = False, tool = 'psfr')
        if ct >8:
            ct = 'all'
        pickle.dump(_data_process , open(run_folder+'fit_material/'+'data_process_idx{0}_{2}_CombPsfsNO_{1}.pkl'.format(idx, ct, filt), 'wb'))
        
    if if_printshow==True:
        fit_run = fit_run_list[idx_counts[0]]
        fit_run.plot_final_qso_fit()
        host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
        AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
        ratio = host_flux/(host_flux+AGN_flux)
        print('count', count)
        print('Chisqs top 2', round(chisqs[idx_counts[0]],2), round(chisqs[idx_counts[1]],2))
        print_s ='idx: '+ item[0]+' '+item[1]+' ratio: ' + str(round(ratio,2)) +" OK?"
        print(fit_files[idx_counts[0]])
        hold = input(print_s)