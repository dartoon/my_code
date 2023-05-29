#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 10:19:04 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
from galight.data_process import DataProcess
from galight.tools.astro_tools import plt_fits
from galight.tools.measure_tools import measure_bkg
import pickle

run_folder = 'stage3_all/' #!!!
filt = 'F356W'
idx = 1

dp_files = glob.glob(run_folder+'fit_material/data_process_idx{1}_{0}_FOVpsf*.pkl'.format(filt,idx))
dp_files.sort()

import sys
sys.path.insert(0,'../model_z6_data_id0/')

from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

#%%

rerun = True

if rerun == True:
    # for i in [0]:
    for i in range(len(dp_files)):
        file = dp_files[i]
        print(file)
        filename = dp_files[i].replace('data_process', 'fit_run_fixn1_').format(i)
        filename = filename.replace('fit_material', 'fit_material_super2').format(i)
        idx = file.split('idx')[1].split('_')[0]
        target_id = target_id
        _data_process = pickle.load(open(file,'rb'))
        psf = _data_process.PSF_list[-1]
        psf[psf<0] = 0.
        psf = abs(psf)
        _data_process.PSF_list[-1] = psf
        _data_process.noise_map = np.nan_to_num(_data_process.noise_map, nan=1000)
        
        fit_sepc = FittingSpecify(_data_process)
        fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3)
                                      # ps_pix_center_list = [ps_pos]  ) #, fix_n_list= [[0,4],[1,1]])
        fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
        fit_sepc.kwargs_params['lens_light_model'][4][0]['R_sersic'] = 1.
        fit_sepc.kwargs_params['lens_light_model'][2][0]['n_sersic'] = 1.
        # fit_sepc.kwargs_constraints['linear_solver'] = False
        fit_sepc.kwargs_numerics['point_source_supersampling_factor'] = 2
        
        fit_sepc.plot_fitting_sets()
        fit_run = FittingProcess(fit_sepc, savename = target_id)
        fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
        # fit_run.plot_final_qso_fit(target_ID =target_id)
        filt = _data_process.filt
        pickle.dump(fit_run , open(filename, 'wb'))


    import copy
    files = glob.glob(run_folder+'fit_material/data_process_idx{1}_*{0}*_*FOVpsf*.pkl'.format(filt,idx))
    files.sort()
    collect_info = []
    for i in range(len(files)):
        _file = files[i]
        idx_info = _file.split('idx')[1].split('_')[0]
        filt_info = _file.split('W_')[0].split('_')[-1] + 'W'
        this_info = [idx_info, filt_info]
        if this_info not in collect_info:
            collect_info.append(this_info)
    
    if_printshow = False
    for count in range(len(collect_info)):
        item = collect_info[count]
        fit_run_list = []
        idx, filt= item
        fit_files = glob.glob(run_folder+'fit_material_super2/fit_run_fixn1__idx{0}_{1}_*FOVpsf*.pkl'.format(idx, filt))
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
            pickle.dump(_data_process , open(run_folder+'fit_material_super2/'+'data_process_idx{0}_{2}_CombPsfsNO_{1}.pkl'.format(idx, ct, filt), 'wb'))
    dp_files = glob.glob(run_folder  + 'fit_material_super2/data_process_idx{1}_*{0}*CombPsfs*.pkl'.format(filt,idx) ) 
    dp_files.sort()
    
    for i in range(len(dp_files)):
        file = dp_files[i]
        print(file)
        filename = dp_files[i].replace('data_process', 'fit_run_fixn1_')[:-4].format(i)
        idx = file.split('idx')[1].split('_')[0]
        target_id = target_id
        _data_process = pickle.load(open(file,'rb'))
        psf = _data_process.PSF_list[-1]
        psf[psf<0] = 0.
        psf = abs(psf)
        _data_process.PSF_list[-1] = psf
        _data_process.noise_map = np.nan_to_num(_data_process.noise_map, nan=1000)
        fit_sepc = FittingSpecify(_data_process)
        fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3)
                                      # ps_pix_center_list = [ps_pos]  ) #, fix_n_list= [[0,4],[1,1]])
        fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
        fit_sepc.kwargs_params['lens_light_model'][4][0]['R_sersic'] = 1.
        fit_sepc.kwargs_params['lens_light_model'][2][0]['n_sersic'] = 1.
        # fit_sepc.kwargs_constraints['linear_solver'] = False
        fit_sepc.kwargs_numerics['point_source_supersampling_factor'] = 2
        fit_sepc.plot_fitting_sets()
        fit_run = FittingProcess(fit_sepc, savename = target_id)
        fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
        fit_run.plot_final_qso_fit(target_ID =target_id)
        filt = _data_process.filt
        pickle.dump(fit_run , open(filename, 'wb'))
    
#%%
idx = 1
for top_psf_id in [0]:
    fit_run_list = []
    fit_files = glob.glob(run_folder+'fit_material_/fit_run_fixn1_*idx{0}_{1}_*.pkl'.format(idx, filt))#+\
    fit_files = glob.glob(run_folder+'fit_material_super2/fit_run_fixn1_*idx{0}_{1}_*.pkl'.format(idx, filt))#+\
    fit_files.sort()
    for i in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort()  
    # print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])

    count_n = 5
    Chisq_best = chisqs[sort_Chisq[top_psf_id]]
    Chisq_last= chisqs[sort_Chisq[count_n-1]]
    inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
    weight = np.zeros(len(chisqs))
    for i in sort_Chisq[:count_n]:
        weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))    
    prop_name = 'magnitude'
    # all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
    all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
    weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
    rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
    
    
    fit_run = fit_run_list[sort_Chisq[top_psf_id]]
    fit_run.plot_final_qso_fit(target_ID =target_id)
    print('reduced_Chisq:', fit_run.reduced_Chisq)
#     # result.append([filt, fit_run.fitting_specify_class.zp, weighted_value, rms_value])
    host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
    AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
    ratio = host_flux/(host_flux+AGN_flux)
    print(ratio)
    print(fit_run.final_result_galaxy[0][prop_name])
    # print(fit_run.final_result_galaxy)
    print(prop_name, round(weighted_value,2), '+-', round(rms_value,2))
    
    # for i in range(5):
    #     print(all_values[sort_Chisq[i]])
    # for i in range(5):
    #     host_flux = fit_run_list[sort_Chisq[i]].final_result_galaxy[0]['flux_within_frame']
    #     AGN_flux = fit_run_list[sort_Chisq[i]].final_result_ps[0]['flux_within_frame']
    #     ratio = host_flux/(host_flux+AGN_flux)
    #     print(ratio)
    
    
    
    