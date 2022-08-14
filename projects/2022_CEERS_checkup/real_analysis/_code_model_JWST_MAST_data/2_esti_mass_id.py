#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 13:24:22 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import sys
import glob
import pickle
from fun_smass_estimate import esti_smass

filters =  ['f115w', 'f150w', 'f200w', 'f277w', 'f356w', 'f410m', 'f444w']
#352, 353, 354, 355, 356, 362, 357

lines = open('fit_result/'+'f356w'+"_targets_info.txt","r").read().split('\n')
target_IDs = [lines[i].split("['")[1].split("'")[0] for i in range(len(lines)) if lines[i] != ''] 

# #Should use the use_PSF_list info to rule out the QSO.
folder = 'output/'
target_id = 'aegis_511'
z = 2.748
pkl_files = glob.glob(folder+'*fit_result_PsfLib_'+target_id+'.pkl')
pkl_files.sort()

mags = []
for filename in pkl_files:
    fit_run_list, use_PSF_list = pickle.load(open(filename,'rb'))
    filt = filename.split('output/')[1].split('_')[0]
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    idx_counts = chisqs.argsort()
    if chisqs[idx_counts[0]] < chisqs[idx_counts[1]]/2.5:
        print("The AGN org image is used PSF in", target_id , filt)
        idx_counts = idx_counts[1:]
    fit_run = fit_run_list[idx_counts[0]]
    # fit_run.plot_final_qso_fit(target_ID = target_id + ' '+filt, cmap = None)
    host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
    AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
    ratio = host_flux/(host_flux+AGN_flux)
    # print(filt, round(fit_run.final_result_galaxy[0]['flux_within_frame'],2),
    #       "host ratio",round(ratio,2),
    #       round(fit_run.final_result_galaxy[0]['magnitude'],2) )
    mags.append(fit_run.final_result_galaxy[0]['magnitude'])

#%%
z = float(z)
if np.max(mags)>0:
    esti_smass('12345', mags, z)