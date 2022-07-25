#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 16:13:13 2022

@author: Dartoon
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

lines = open(folder+'f356w'+"_targets_info.txt","r").read().split('\n')
target_IDs = [lines[i].split("['")[1].split("'")[0] for i in range(len(lines)) if lines[i] != ''] 
    # ID = target_IDs[10]
    # print(ID)
for ID in target_IDs:
    filts = ['f150w', 'f356w']
    for filt in filts:
        filename = folder +filt+'_fit_result_PsfLib'+'_'+ID+'.pkl'
        fit_run_list, use_PSF_list = pickle.load(open(filename,'rb'))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        idx_counts = chisqs.argsort()
        if chisqs[idx_counts[0]] < chisqs[idx_counts[1]]/2.5:
            print("The AGN org image is used PSF in",ID , filt)
            idx_counts = idx_counts[1:]
        fit_run = fit_run_list[idx_counts[0]]
        # fit_run.savename = ID +'_' +filt
        fit_run.plot_final_qso_fit(target_ID = ID +'_' +filt, save_plot=False)
        host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
        AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
        ratio = host_flux/(host_flux+AGN_flux)
        print(filt, round(fit_run.final_result_galaxy[0]['flux_within_frame'],2),
              "host ratio",round(ratio,2),
              round(fit_run.final_result_galaxy[0]['magnitude'],2) )

