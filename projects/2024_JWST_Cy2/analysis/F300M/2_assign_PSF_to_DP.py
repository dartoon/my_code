#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 11:47:01 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import glob
# files = glob.glob('material/data_process+apertures_*.pkl')

idx = 1 #!!!
import sys
sys.path.insert(0, '../../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
run_folder = 'stage3_all/' #!!!

# for filt in ['F356W', 'F150W'] :
for filt in ['F300M'] :
    data_process = pickle.load(open(run_folder+'material/data_process_idx{0}_{1}.pkl'.format(idx, filt),'rb'))
    # print(data_process.target_stamp[40,40])
    filt = data_process.filt
    
    PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
    PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
    for i in range(len(PSF_list_clean)):
        psf = PSF_list_clean[i]
        if filt == 'F150W':
            psf = psf[40:-40, 40:-40]
        data_process.PSF_list = [psf]
        pickle.dump(data_process , open(run_folder+'fit_material/'+'data_process_idx{2}_{0}_FOVpsf{1}.pkl'.format(filt, i, idx), 'wb'))
        
    # PSF_lib_files = glob.glob('../../2022_CEERS_checkup/real_analysis/model_JWST_stage3_Masafusa/material/*'+filt[:-1]+'*_PSF_Library.pkl')[0]
    # PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
    # for i in range(len(PSF_list_clean)):
    #     psf = PSF_list_clean[i]
    #     data_process.PSF_list = [psf]
    #     pickle.dump(data_process , open(run_folder+'fit_material/'+'data_process_idx{2}_{0}_CEERSpsf{1}.pkl'.format(filt, i, idx), 'wb'))
        
