#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 11:43:09 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle
run_folder = '../model_z6_*/stage3_all/' #!!!
import glob
# for filt in ['F356W', 'F150W'] :
    
#%% SNR of the QSO image:    
for idx in [0,1]:
    for filt in ['F356W', 'F150W'] :
        files = glob.glob(run_folder+'material/data_process_idx{0}_{1}.pkl'.format(idx, filt))
        print(files)
        file = files[0]
        data_process = pickle.load(open(file,'rb'))
        print(np.sqrt(np.sum((data_process.target_stamp[20:-20,20:-20]/data_process.noise_map[20:-20,20:-20])**2)))
#%% SNR of the QSO image:    
for idx in [0,1]:
    for filt in ['F356W', 'F150W'] :
            top_psf_id =0
            fit_run_list = []
            # idx = idx_info
            if idx == 0:
                fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_FOV*.pkl'.format(idx, filt))#+\
            elif idx == 1:
                fit_files = glob.glob(run_folder+'*fit_material*/fit_run_fixn1__idx{0}_{1}_FOV*.pkl'.format(idx, filt))#+\
            fit_files.sort()
            for i in range(len(fit_files)):
                fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
            chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
            sort_Chisq = chisqs.argsort()  
            fit_run = fit_run_list[sort_Chisq[top_psf_id]]
            print(idx, filt)
            # print(np.sqrt(np.sum((fit_run.flux_2d_out['data-Point Source']/
            #                       fit_run.fitting_specify_class.data_process_class.noise_map)**2)))
            print(np.sqrt(np.sum((fit_run.image_host_list[0]/
                                  fit_run.fitting_specify_class.data_process_class.noise_map)**2)))
