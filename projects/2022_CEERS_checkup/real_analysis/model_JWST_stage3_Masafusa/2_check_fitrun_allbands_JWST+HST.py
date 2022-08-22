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

remove_id = [24, 55]

# for idx in range(11,70):
# for idx in [0, 2, 8, 28, 35, 51]:
for idx in [31, 32, 56]:
    if idx in remove_id:
        continue
    else:
        files = glob.glob('../*/fit_material/data_process_idx{0}_*_psf*.pkl'.format(idx))
        files.sort()
        
        collect_info = []
        for i in range(len(files)):
            _file = files[i]
            idx_info = _file.split('idx')[1].split('_')[0]
            filt_info = _file.split('_psf')[0].split('_')[-1]
            this_info = [idx_info, filt_info]
            if this_info not in collect_info:
                collect_info.append(this_info)
        
        f = open("../model_JWST_stage3_Masafusa/target_idx_info.txt","r")
        string = f.read()
        lines = string.split('\n')   # Split in to \n
        target_id = [lines[i].split(' ')[1] for i in range(len(lines)) if lines[i].split(' ')[0] == str(idx)][0]
        
        for count in range(len(collect_info)):
            item = collect_info[count]
            fit_run_list = []
            idx, filt= item
            fit_files = glob.glob('../model_JWST_stage3_Masafusa/fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))+\
                        glob.glob('../model_HST_material/fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))
            fit_files.sort()
            warn_strs = ['F115W_psf6', 'F150W_psf7', 'F277W_psf2']
            for warn_str in warn_strs:
                fit_files = [fit_files[i] for i in range(len(fit_files)) if warn_str not in fit_files[i]]
            
            for i in range(len(fit_files)):
                fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
            chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
            sort_Chisq = chisqs.argsort()  
            print('idx', idx, filt, "Total PSF NO.", len(sort_Chisq))
            
            if 'HST_material' in fit_files[0]:
                weight = np.zeros(len(chisqs))
                weight[sort_Chisq[0]] = 1
            elif 'JWST_stage3' in fit_files[0]:
                count_n = 5
                Chisq_best = chisqs[sort_Chisq[0]]
                Chisq_last= chisqs[sort_Chisq[count_n-1]]
                inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
                weight = np.zeros(len(chisqs))
                for i in sort_Chisq[:count_n]:
                    weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
            fit_run = fit_run_list[sort_Chisq[0]]
            fit_run.plot_final_qso_fit(target_ID = target_id+'-'+filt)
            
            prop_name = 'n_sersic'
            all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
            weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
            rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
            
            host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
            AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
            ratio = host_flux/(host_flux+AGN_flux)
            print(prop_name, round(weighted_value,2), '+-', round(rms_value,2))
            print('Chisqs top 2', round(chisqs[sort_Chisq[0]],2), round(chisqs[sort_Chisq[1]],2))
            print_s ='idx: '+ item[0]+' '+item[1]+' ratio: ' + str(round(ratio,2)) +" OK? \n\n\n"
            print(fit_files[sort_Chisq[0]])
            print(print_s)
            # hold = input(print_s)
        hold = input("idx {0} above, OK?\n\n".format(idx))