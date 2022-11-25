#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 15:03:04 2022

@author: Dartoon
"""


import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle
import copy


run_folder = 'stage3_all/' #!!!
filt = 'F356W'
files = glob.glob(run_folder+'fit_material/data_process_idx0_*{0}*_*FOVpsf*.pkl'.format(filt))
files.sort()
collect_info = []
for i in range(len(files)):
    _file = files[i]
    idx_info = _file.split('idx')[1].split('_')[0]
    filt_info = _file.split('W_')[0].split('_')[-1] + 'W'
    this_info = [idx_info, filt_info]
    if this_info not in collect_info:
        collect_info.append(this_info)

# - [ ] F115W_psf6 is QSO (idx 2). [136, 137, 139]
# - [ ] F150W_psf7 is QSO (idx 2).
# - [ ] F277W_psf2 is QSO (idx 2).

print("After remove candidates")
PSF_lib_files = glob.glob('stage3_all/'+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(0))[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))

#%%
if_printshow = False
for count in range(len(collect_info)):
    item = collect_info[count]
    fit_run_list = []
    idx, filt= item
    fit_files = glob.glob(run_folder+'fit_material/fit_run_idx{0}_{1}_*FOVpsf*.pkl'.format(idx, filt))
    fit_files.sort()
    
    for i in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    idx_counts = chisqs.argsort()  
    if len(idx_counts)<8:
        print(idx, filt, len(idx_counts))
    # print("work on", count, 'idx', idx, filt, "Total PSF NO.", len(idx_counts))
    for i in range(5):
        print(PSF_RA_DEC_list[idx_counts[i]])