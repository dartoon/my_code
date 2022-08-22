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

# # for file in glob.glob("fit_material/fit_run_idx10_F150W_*")[:]:
# for file in ['fit_material/fit_run_idx10_F115W_CombPsfsNO_8_10.pkl']:
#     fit_run = pickle.load(open(file,'rb'))
#     # fit_run.plot_final_qso_fit()
#     print(fit_run.final_result_ps[0]['point_amp'])
#     print(fit_run.fitting_seq.likelihoodModule.log_likelihood(verbose=True, kwargs_return=fit_run.kwargs_result))

#%%
# f = open("../model_JWST_stage3_Masafusa/target_idx_info.txt","r")
# string = f.read()
# lines = string.split('\n')   # Split in to \n

# dp_files = glob.glob('fit_material/data_process_idx*.pkl')
# dp_files.sort()
# for i in range(len(dp_files)):
files = glob.glob('fit_material/'+'fit_run_idx*.pkl')
files.sort()
check_list = []
for i in range(0,len(files)):
    file = files[i]
    fit_run = pickle.load(open(file,'rb'))
    host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
    AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
    ratio = host_flux/(host_flux+AGN_flux)
    idx_info = file.split('idx')[1].split('_')[0]
    if ratio > 1:
        check_list.append([idx_info, file, ratio])
    
# #%%Reload fit_run and perform fitting to avoid negative PS.
# from galight.fitting_process import FittingProcess
# from galight.fitting_specify import FittingSpecify
# import shutil
# from galight.tools.measure_tools import mask_obj
# # for i, item in enumerate(check_list[:100]):
# for i, item in enumerate(check_list[100:200]):
# # for i, item in enumerate(check_list[200:]):
#     _, filename,_ = item
#     fit_run = pickle.load(open(filename,'rb'))
#     _data_process = fit_run.fitting_specify_class.data_process_class
#     _data_process.PSF_list[-1] = abs(_data_process.PSF_list[-1])
#     fit_sepc = FittingSpecify(_data_process)
#     target_stamp = _data_process.target_stamp 
#     mask = mask_obj(target_stamp, _data_process.apertures[:1], if_plot=False)
    
#     ps_pos = np.where(target_stamp == np.max(target_stamp * (1-mask[0])))
#     ps_pos = (ps_pos[0][0] - _data_process.radius, ps_pos[1][0] - _data_process.radius)
#     ps_pos = [ps_pos[1], ps_pos[0]]
    
#     fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3, ps_pix_center_list = [ps_pos]) #, fix_n_list= [[0,4],[1,1]])
#     fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
#     fit_sepc.kwargs_constraints['linear_solver'] = False
#     fit_sepc.build_fitting_seq()
#     fit_sepc.kwargs_params['lens_light_model'][0] = fit_run.kwargs_result['kwargs_lens_light']
#     fit_sepc.plot_fitting_sets()
#     fit_run = FittingProcess(fit_sepc, savename = None)
#     fit_run.run(algorithm_list = ['PSO','PSO'], fitting_level=['norm','deep'])
#     bk_filename = filename.replace('fit_material', 'backup_files')
#     shutil.move(filename, bk_filename)
#     pickle.dump(fit_run , open(filename, 'wb'))
#     fit_run.plot_final_qso_fit()
#     print(i, 'finished')
    

# #%%
# files = glob.glob('fit_material/data_process_idx*_psf*.pkl')
# files.sort()
# collect_info = []
# for i in range(len(files)):
#     _file = files[i]
#     idx_info = _file.split('idx')[1].split('_')[0]
#     filt_info = _file.split('_psf')[0].split('_')[-1]
#     this_info = [idx_info, filt_info]
#     if this_info not in collect_info:
#         collect_info.append(this_info)

# if_printshow = True
# for count in range(len(collect_info)):
#     item = collect_info[count]
#     fit_run_list = []
#     idx, filt= item
#     fit_files = glob.glob('fit_material/fit_run_idx{0}_{1}_C*.pkl'.format(idx, filt))
#     fit_files.sort()
#     warn_strs = ['F115W_psf6', 'F150W_psf7', 'F277W_psf2']
#     for warn_str in warn_strs:
#         fit_files = [fit_files[i] for i in range(len(fit_files)) if warn_str not in fit_files[i]]
    
#     for i in range(len(fit_files)):
#         fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
#     chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
#     idx_counts = chisqs.argsort()  
#     if len(idx_counts)<8:
#         print(idx, filt, len(idx_counts))
#     print("work on", count, 'idx', idx, filt, "Total PSF NO.", len(idx_counts))
        
#     if if_printshow==True:
#         fit_run = fit_run_list[idx_counts[0]]
#         fit_run.plot_final_qso_fit()
#         host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
#         AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
#         ratio = host_flux/(host_flux+AGN_flux)
#         print('count', count)
#         print('Chisqs top 2', round(chisqs[idx_counts[0]],2), round(chisqs[idx_counts[1]],2))
#         print_s ='idx: '+ item[0]+' '+item[1]+' ratio: ' + str(round(ratio,2)) +" OK?"
#         print(fit_files[idx_counts[0]])
#         hold = input(print_s)