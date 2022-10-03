#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 13:23:13 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from functions_for_result import esti_smass, load_prop, load_info, name_list
import glob


#%%
import glob, pickle
# idx = 1
# filt = 'F444W'
def cal_mag_err(idx, filt):
    files = glob.glob('../*/fit_material/data_process_idx{0}_*_psf*.pkl'.format(idx))
    files.sort()
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
    # _, z = load_info(idx)
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
    fit_run = fit_run_list[sort_Chisq[0]]
    count_n = 5
    Chisq_best = chisqs[sort_Chisq[0]]
    Chisq_last= chisqs[sort_Chisq[count_n-1]]
    inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
    weight = np.zeros(len(chisqs))
    for i in sort_Chisq[:count_n]:
        weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
        
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
    
    prop_name = 'magnitude'
    all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
    weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
    rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
    if idx in ['35', '0'] and filt == 'F444W':  #!!!]
        rms_value = 0.4
    if idx in ['2'] and filt == 'F410M':  #!!!]
        rms_value = 0.2
    return rms_value


# ID, mags, z = 'idx0', 
# 1,2,0,51,35
# folder = '202209' #0.2 mag no HST.
# folder = '20220901_' #Not HST
folder = '20221003' #HST upper limit
rerun = True
# idx = [1,2,0,51,35]
# for idx in [35,0,2,51,1]:
for idx in [1]:
    target_id, z = load_info(idx)
    if rerun == True:
        if idx == 1:
            root_folder = '../*/*'  #Include HST
        else:
            root_folder = './*'  #Include HST
        mag_result = load_prop(idx, root_folder = root_folder)
        if idx ==1 and 'F814W' in mag_result.keys():
            mag_result['F814W'] += -2.5*np.log10(5) 
            del mag_result['F606W']
            del mag_result['F125W']
            del mag_result['F160W']
        import time
        t1 = time.time()
        print('Run estimate')
        band_as_upper = []
        if idx ==35:
            band_as_upper = ['F115W']
        if idx ==2:
            band_as_upper = ['F115W']
            mag_result['F115W'] += -2.5*np.log10(5)
        
        mag_err = []
        for filt in mag_result.keys():
            err = cal_mag_err(str(idx), filt)
            if err<0.15:
                err = 0.15
            mag_err.append( err )
        esti_smass(ID = folder+str(idx), mags_dict = mag_result, z = z, flag = 0, 
                    if_run_gsf=True, band_as_upper = band_as_upper, mag_err=mag_err)
        t2 = time.time()
