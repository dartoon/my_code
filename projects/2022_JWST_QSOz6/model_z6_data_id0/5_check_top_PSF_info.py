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


filt = 'F150W'
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
#%%
PSF_lib_files = glob.glob('stage3_all/'+'material/*'+filt[:-1]+'*_PSF_Library_idx?.pkl')[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
#%%
if_printshow = False
item = collect_info[0]
fit_run_list = []
idx, filt= item
fit_files = glob.glob(run_folder+'fit_material/fit_run_idx{0}_{1}_*FOVpsf*.pkl'.format(idx, filt))
fit_files.sort()
for i in range(len(fit_files)):
    fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
idx_counts = chisqs.argsort()  
ct = 3
PSF_RA_DEC_info = [PSF_RA_DEC_list[i] for i in idx_counts[:ct]] 
zp = 27.980780691581828
PSF_F356W_mag = [-2.5*np.log10(np.sum(PSF_list_clean[i]))+zp for i in idx_counts[:ct]] 


#%%
PSF_lib_files = glob.glob('stage3_all/'+'material/*'+'150'+'*_PSF_Library_idx?.pkl')[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
PSF_RA_DEC_list = np.array(PSF_RA_DEC_list)
PSF_F150W_mag = [None]*ct
    
zp = 28.03341727868797
for i in range(len(PSF_RA_DEC_info)):
    dis = np.sqrt(np.sum((PSF_RA_DEC_info[i] - PSF_RA_DEC_list)**2,axis = 1))*3600
    if np.min(dis) < 5*0.03:
        psf_id = np.where(dis ==dis.min())[0][0]
        PSF_F150W_mag[i] =  -2.5*np.log10(np.sum(PSF_list_clean[psf_id][40:-40, 40:-40])) + zp
        
print(PSF_RA_DEC_info, PSF_F356W_mag, PSF_F150W_mag)