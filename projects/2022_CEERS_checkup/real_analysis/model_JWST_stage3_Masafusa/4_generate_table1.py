#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 20:30:48 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import os
import glob
from galight.tools.asymmetry_tools import CAS

remove_id = [24, 55]

#!!! 0 141943.58+525431.3 214.93160692630738 52.90871362182564 3.442 -99
#!!! 1 142005.59+530036.7 215.0233219305857 53.01020458328741 1.646 -99
#!!! 2 142008.61+530004.0 215.03590035000232 53.00111935271983 2.588 -99

#### 8 aegis_585 214.93166 52.908708 3.435 3.364  #Same as idx 0
#### 26 aegis_742 215.02339 53.010207 1.644 1.74  # Same as SDSS idx 1
#### 28 aegis_463 214.7768 52.825876 2.274 2.24  #A type 2
#!!! 35 aegis_482 214.75522 52.836795 3.465 3.345
#!!! 51 aegis_477 214.87073 52.833117 2.317 2.136

#### 29 aegis_465 214.74245 52.826187 -99.0 3.019 #Photoz>3
#### 53 aegis_495 214.87124 52.845067 -99.0 3.422 #Photoz>3

# for idx in range(1,2):
# for idx in [31, 32, 56]:   #JWST PS position corrected
# for idx in [0, 2, 8, 28, 35, 51]:  #z_spec > 2
# for idx in [29, 53]:  #z_spec > 2


# for idx in [0, 1, 2, 35, 51]:  #z_spec > 1.6, QSO 
result_list = []
target_ID_list = []
filt = 'F444W'
cal_CAS = True
if_plot = False
# for idx in [1,2,0,51,35]:  #z_spec > 1.6
for idx in [35,0,2,51,1]:  #z_spec > 1.6
# for idx in [1]:  #z_spec > 1.6
    result = []
    files = glob.glob('../*/*fit_material/data_process_idx{0}_*_psf*.pkl'.format(idx))
    files.sort()
    # file_ACS = glob.glob('../model_HST_material/material/data_process+apertures_{0}_ACS.pkl'.format(idx))
    # file_WFC= glob.glob('../model_HST_material/material/data_process+apertures_{0}_IR.pkl'.format(idx))
    # file_JWST = glob.glob('../model_JWST_stage3_Masafusa/material/data_process+apertures_{0}.pkl'.format(idx))
    
    _collect_info = []
    for i in range(len(files)):
        _file = files[i]
        idx_info = _file.split('idx')[1].split('_')[0]
        filt_info = _file.split('_psf')[0].split('_')[-1]
        this_info = [idx_info, filt_info]
        if this_info not in _collect_info:
            _collect_info.append(this_info)
    
    filters = [_collect_info[i][1] for i in range(len(_collect_info))]
    if 'F814W' in filters:
        filters = ['F814W'] + [filters[i] for i in range(len(filters)) if filters[i]!='F814W']
    if 'F606W' in filters:
        filters = ['F606W'] + [filters[i] for i in range(len(filters)) if filters[i]!='F606W']
    
        
    f = open("../model_JWST_stage3_Masafusa/target_idx_info.txt","r")
    string = f.read()
    lines = string.split('\n')   # Split in to \n
    target_id = [lines[i].split(' ')[1] for i in range(len(lines)) if lines[i].split(' ')[0] == str(idx)][0]
    
    fit_run_list = []
    idx = idx_info
    # idx, filt= item
    fit_files = glob.glob('../model_JWST_stage3_Masafusa/*fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))+\
                glob.glob('../model_HST_material/*fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))
    fit_files.sort()
    warn_strs = ['F115W_psf6', 'F150W_psf7', 'F277W_psf2']
    for warn_str in warn_strs:
        fit_files = [fit_files[i] for i in range(len(fit_files)) if warn_str not in fit_files[i]]
    
    for i in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort()  
    print('idx', idx, filt, "Total PSF NO.", len(sort_Chisq))
    target_ID_list.append(target_id)
    if fit_files == []:
        result_list.append([-99]*8)
        continue
    # if 'HST_material' in fit_files[0]:
    #     weight = np.zeros(len(chisqs))
    #     weight[sort_Chisq[0]] = 1
    # elif 'JWST_stage3' in fit_files[0]:
    #     count_n = 5
    #     Chisq_best = chisqs[sort_Chisq[0]]
    #     Chisq_last= chisqs[sort_Chisq[count_n-1]]
    #     inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
    #     weight = np.zeros(len(chisqs))
    #     for i in sort_Chisq[:count_n]:
    #         weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
    # prop_name = 'R_sersic'
    # # all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
    # all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
    # weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
    # rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
    
    
    fit_run = fit_run_list[sort_Chisq[0]]
    if idx == '35' and filt == 'F444W':  #!!!
        fit_run = fit_run_list[sort_Chisq[3]]
    if idx == '0' and filt == 'F444W':  #!!!
        fit_run = fit_run_list[sort_Chisq[1]]
    if if_plot == True:
        fit_run.plot_final_qso_fit(target_ID = target_id+'$-$'+filt)
    host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
    AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
    ratio = host_flux/(host_flux+AGN_flux)
    if filt == 'F444W':
        if fit_run.fitting_specify_class.data_process_class.target_pos[0] < 5000: #module A
            correct = 0.44157708 / 0.343
            mag_correct = +2.5*np.log10(correct)
        if fit_run.fitting_specify_class.data_process_class.target_pos[0] > 5000: #module B
            correct = 0.3899884 / 0.335
            mag_correct = +2.5*np.log10(correct)
        fit_run.final_result_galaxy[0]['magnitude'] = fit_run.final_result_galaxy[0]['magnitude'] + mag_correct
    if filt == 'F410M':
        if fit_run.fitting_specify_class.data_process_class.target_pos[0] < 5000:
            correct = 0.9355298 / 0.832
            mag_correct = +2.5*np.log10(correct)
        if fit_run.fitting_specify_class.data_process_class.target_pos[0] > 5000:
            correct = 0.9272488 / 0.811
            mag_correct = +2.5*np.log10(correct)
        fit_run.final_result_galaxy[0]['magnitude'] = fit_run.final_result_galaxy[0]['magnitude'] + mag_correct
    
    if cal_CAS == True:
    
        CAS_class = CAS(fit_run, seg_cal_reg = 'or', obj_id=0, rm_ps=True, rm_obj=True)
        cas = CAS_class.cal_CAS(mask_type='segm',extend=1, if_plot=False, radius=35)
        # print('asymmetry:', cas[0])
        # print('smoothness (by abs diff and pos diff):', cas[1])
        # print('concentration:', cas[2])
        # print('Gini:', cas[3])
    else:
        cas = 0, [0,0], 0
    result.append([fit_run.reduced_Chisq, ratio * 100, fit_run.final_result_galaxy[0]['R_sersic'],
                   fit_run.final_result_galaxy[0]['n_sersic'], fit_run.final_result_galaxy[0]['magnitude'],
                   cas[2], cas[1][0], cas[0]])
    # print(prop_name, round(weighted_value,2), '+-', round(rms_value,2))
    # print('Chisqs top 2', round(chisqs[sort_Chisq[0]],2), round(chisqs[sort_Chisq[1]],2))
    # print_s =filt +' ratio: ' + str(round(ratio,2)) + "\n\n\n"
    print(fit_files[sort_Chisq[0]])
    # print(print_s)
    # hold = input(print_s)
    result_list.append(result[0])

#%%print:
prop_list = ['$\chi ^2$ (Reduced)', 'host-total flux ratio', 'Reff ($\\arcsec$)',
             '\sersic\ $n$', 'host mag', 'concentration', 'asymmetry','smoothness']
digt_list = [2, 1, 3, 1, 2, 2,2,2]
for i in range(len(result_list[0])):
    print_material = [result_list[j][i] for j in range(len(result_list)) ]
    print_s = '& '
    for j in range(len(print_material)):
        if j != 0:
            print_s = print_s + ' & '
        if i != 1:
            print_s = print_s + str(round(print_material[j],digt_list[i]))
        elif i == 1:
            print_s = print_s + str(round(print_material[j],digt_list[i])) + '\% '
    print_s = print_s + '\\\\'
    print(prop_list[i], print_s)
