#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 11:25:37 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle
import copy

files = glob.glob('fit_material/data_process_idx*.pkl')
files.sort()
collect_info = []
for i in range(len(files)):
    _file = files[i]
    idx_info = _file.split('idx')[1].split('_')[0]
    filt_info = _file.split('_psf')[0].split('_')[-1]
    this_info = [idx_info, filt_info]
    if this_info not in collect_info:
        collect_info.append(this_info)

#%%
for count in range(len(collect_info)):
# for count in [136]:
    item = collect_info[count]
    fit_run_list = []
    idx, filt= item
    fit_files = glob.glob('fit_material/fit_run_idx{0}_{1}_psf*.pkl'.format(idx, filt))
    fit_files.sort()
    for i in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    idx_counts = chisqs.argsort()  
    fit_run = fit_run_list[idx_counts[0]]
    # fit_run.plot_final_qso_fit()
    host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
    AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
    ratio = host_flux/(host_flux+AGN_flux)
    print(fit_files[idx_counts[0]])
    print(round(chisqs[idx_counts[0]],2), round(chisqs[idx_counts[1]],2))
    print(count)
    print_s = item[0]+' '+item[1]+ ' Chisq: '+ str(chisqs[idx_counts[0]]) +' ratio: ' + str(round(ratio,2)) +" OK?"
    hold = input(print_s)
    # print(print_s)
    # if len(idx_counts)<8:
    #     print(idx, filt, len(idx_counts))
    # print
    # for ct in [3, 5, 8, 'all']:
    #     _data_process_list = [fit_run_list[i].fitting_specify_class.data_process_class for i in  idx_counts[:5]] 
        # if ct != 'all':
        #     PSF_list_for_comb = [PSF_list_clean[i] for i in idx_counts[:ct]]
        # elif ct == 'all':
        #     PSF_list_for_comb = [PSF_list_clean[i] for i in idx_counts]
        # _data_process.PSF_list  = copy.deepcopy(PSF_list_for_comb)
        # _data_process.stack_PSF(if_plot = True, tool = 'psfr')