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

#!!! 0 141943.58+525431.3 214.93160692630738 52.90871362182564 3.442 -99
#!!! 1 142005.59+530036.7 215.0233219305857 53.01020458328741 1.646 -99
#!!! 2 142008.61+530004.0 215.03590035000232 53.00111935271983 2.588 -99

#### 8 aegis_585 214.93166 52.908708 3.435 3.364  #Same as idx 0
#### 26 aegis_742 215.02339 53.010207 1.644 1.74  # Same as SDSS idx 1
#### 28 aegis_463 214.7768 52.825876 2.274 2.24  #A type 2
#!!! 35 aegis_482 214.75522 52.836795 3.465 3.345
#!!! 51 aegis_477 214.87073 52.833117 2.317 2.136

#### 29 aegis_465 214.74245 52.826187 -99.0 3.019 #Photoz>3
#### 53 aegis_495 214.87124 52.845067 -99.0 3.422 #Photoz>3  #!!! In CEERS
#### 57 aegis_511 214.89561 52.856516  !!! In CEERS
#### 39 aegis_525 214.85389 52.861419  !!! In CEERS
#### 41 aegis_532 214.8506 52.866459  !!! In CEERS

# for idx in range(1,2):
# for idx in [31, 32, 56]:   #JWST PS position corrected
# for idx in [0, 2, 8, 28, 35, 51]:  #z_spec > 2
# for idx in [29, 53]:  #z_spec > 2
from functions_for_result import name_list
result = []
for idx in range(10, 21):  #CEERS
    if idx in remove_id:
        continue
    else:
        # files = glob.glob('../*/*fit_material/data_process_idx{0}_*_psf*.pkl'.format(idx))
        files = glob.glob('./*fit_material/data_process_idx{0}_*_psf*.pkl'.format(idx))
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
            filters = ['F814W'] + [filters[i] for i in range(len(filters)) if filters[i]!='F814W' ]
        if 'F606W' in filters:
            filters = ['F606W'] + [filters[i] for i in range(len(filters)) if filters[i]!='F606W' ]
        
            
        f = open("../model_JWST_stage3_Masafusa/target_idx_info.txt","r")
        string = f.read()
        lines = string.split('\n')   # Split in to \n
        target_id = [lines[i].split(' ')[1] for i in range(len(lines)) if lines[i].split(' ')[0] == str(idx)][0]
        
        f = open("./material/target_info.txt","r")
        string_1 = f.read()
        lines_ = string_1.split('\n')   # Split in to \n
        spec_z = [lines_[i].split(' ')[3] for i in range(len(lines_)) if lines_[i].split(' ')[0] == target_id][0]
        photo_z = [lines_[i].split(' ')[4] for i in range(len(lines_)) if lines_[i].split(' ')[0] == target_id][0]
        
        if float(spec_z) >0:
            z_str = 'spec_z=' + spec_z
        else:
            z_str = 'photo=' + photo_z
        
        filters = [filters[0], filters[-1]]
        
        for count in range(len(filters)):
            fit_run_list = []
            idx = idx_info
            filt = filters[count]
            # idx, filt= item
            fit_files = glob.glob('./*fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
                        # glob.glob('../model_HST_material/fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))
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
            else:
                count_n = 5
                Chisq_best = chisqs[sort_Chisq[0]]
                Chisq_last= chisqs[sort_Chisq[count_n-1]]
                inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
                weight = np.zeros(len(chisqs))
                for i in sort_Chisq[:count_n]:
                    weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
            fit_run = fit_run_list[sort_Chisq[0]]
            if idx == '35' and filt == 'F444W':  #!!!
                fit_run = fit_run_list[sort_Chisq[3]]
            if idx == '0' and filt == 'F444W':  #!!!
                fit_run = fit_run_list[sort_Chisq[1]]
            
            fit_run.savename = 'fit_collection/{0}_{1}_{2}.pdf'.format(idx, target_id, filt)
            fit_run.plot_final_qso_fit(target_ID = target_id+'$-$'+filt + z_str, save_plot = True)
            
            # name_list = {1:'SDSS1420+5300A', 2: 'SDSS1420+5300B', 0: 'SDSS1419+5w254', 
            #               51: 'aegis_477', 35: 'aegis_482'}

            # target_id = name_list[int(idx)]
            # fit_run.savename = 'outcomes/'+'ID' + idx + '_' + filt
            # fit_run.plot_final_qso_fit(target_ID = target_id+'  '+filt, save_plot= True)
            
            prop_name = 'magnitude'
            # all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
            all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
            weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
            rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
            
            result.append([filt, fit_run.fitting_specify_class.zp, weighted_value, rms_value])
            
            host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
            AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
            ratio = host_flux/(host_flux+AGN_flux)
            print(fit_run.final_result_galaxy)
            print(prop_name, round(weighted_value,2), '+-', round(rms_value,2))
            print('Chisqs top 2', round(chisqs[sort_Chisq[0]],2), round(chisqs[sort_Chisq[1]],2))
            print_s =filt +' ratio: ' + str(round(ratio,2)) + "\n\n\n"
            print(fit_files[sort_Chisq[0]])
            print(print_s)
            # hold = input(print_s)
        # hold = input("idx {0} above, OK?\n\n".format(idx))