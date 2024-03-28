#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 21:30:09 2024

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle

idx = 5
filt = 'F356W'
run_folder = '../material/fit_result/'
psf_sp = 2
# psf_idx = 1
check_files = []
for psf_idx in range(10):
    load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx{2}*_psfsf{3}.pkl'.format(filt, idx, psf_idx, psf_sp))
    load_files.sort()
    chisqs_idx = []
    for file in load_files:
        fit_run = pickle.load(open(file ,'rb'))
        chisqs_idx.append(fit_run.reduced_Chisq)
        if psf_idx == idx:
            print(fit_run.final_result_galaxy[0]['flux_within_frame'], fit_run.reduced_Chisq)
    
    # print(psf_idx, min(chisqs_idx), np.where(min(chisqs_idx) == chisqs_idx)[0][0])
    psf_i = np.where(min(chisqs_idx) == chisqs_idx)[0][0]
    check_files.append(load_files[psf_i] ) 

#%%
print('Check Overall')
for i, file in enumerate(check_files[:10]):
    fit_run = pickle.load(open(file ,'rb'))
    fit_run.plot_final_qso_fit()
    if idx == i:
        print(fit_run.final_result_galaxy[0]['flux_within_frame'], fit_run.reduced_Chisq,'The one')
    else:
        print(fit_run.final_result_galaxy[0]['flux_within_frame'], fit_run.reduced_Chisq)
        