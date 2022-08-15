#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 14:46:03 2022

@author: Dartoon

Load pickles
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle

folder = 'fit_result/'

# ID = '142008.61+530004.0'
# ID = 'aegis_630'
# filt = 'f356w'
# filt = 'f150w'

#Should use the use_PSF_list info to rule out the QSO.

# lines = open(folder+'f356w'+"_targets_info.txt","r").read().split('\n')
# target_IDs = [lines[i].split("['")[1].split("'")[0] for i in range(len(lines)) if lines[i] != ''] 
    # ID = target_IDs[10]
    # print(ID)
# for ID in target_IDs:
for ID in ['142008.61+530004.0']:
    filts = ['F277W','F356W', 'F410M', 'F444W'][:]
    for filt in filts:
        # filename = folder +filt+'_top8PSFcomb_PsfLib'+'_'+ID+'.pkl'
        filename = folder +filt+'_PSFcomb_OrgNoiseMap_PsfLib'+'_'+ID+'.pkl'
        fit_run_list, use_PSF_list = pickle.load(open(filename,'rb'))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        idx_counts = chisqs.argsort()
        if chisqs[idx_counts[0]] < chisqs[idx_counts[1]]/2.5:
            print("The AGN org image is used PSF in",ID , filt)
            idx_counts = idx_counts[1:]
        idx = idx_counts[1]
        fit_run = fit_run_list[idx]
        # fit_run.savename = ID +'_' +filt
        fit_run.plot_final_qso_fit(target_ID = ID +'_' +filt, save_plot=False)
        print(fit_run.reduced_Chisq)
        host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
        AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
        ratio = host_flux/(host_flux+AGN_flux)
        print(use_PSF_list[idx])
        print(filt, '\n', round(fit_run.final_result_galaxy[0]['flux_within_frame'],2),
              "host ratio",round(ratio,2),
              "mag",round(fit_run.final_result_galaxy[0]['magnitude'],2),
              "n_sersic",round(fit_run.final_result_galaxy[0]['n_sersic'],2),
              "R_sersic",round(fit_run.final_result_galaxy[0]['R_sersic'],2),
              )

        print([[round(fit_run_list[i].final_result_galaxy[0]['R_sersic'],2), 
               round(fit_run_list[i].reduced_Chisq, 2)] for i in 
               range(len(fit_run_list))])