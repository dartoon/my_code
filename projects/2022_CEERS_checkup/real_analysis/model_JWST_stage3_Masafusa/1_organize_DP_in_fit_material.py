#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:10:16 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle
# pickle.dump([[data_process_list_l, com_aper_l], [data_process_list_s, com_aper_s] ], open('material/'+'data_process+apertures_{0}.pkl'.format(idx), 'wb'))

# for idx in files
import glob
# files = glob.glob('material/data_process+apertures_*.pkl')

f = open("material/target_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
result_folder = 'fit_result/'
lines = lines[1:]

remove_id = [24, 55]
for idx in range(59):    
    if idx in remove_id:
        continue
    line = lines[idx]
    target_id, RA, Dec, spec_z, photo_z = line.split(' ')
    print(idx, target_id, RA, Dec)
    RA, Dec, spec_z, photo_z = float(RA), float(Dec), float(spec_z), float(photo_z)
    JWST_data_process_list = pickle.load(open('material/data_process+apertures_{0}.pkl'.format(idx),'rb'))
    [data_process_list_l, com_aper_l], [data_process_list_s, com_aper_s]  = JWST_data_process_list
    
    for data_process in data_process_list_l:
        data_process.apertures = com_aper_l
        filt = data_process.filt
        PSF_lib_files = glob.glob('material/*'+filt[:-1]+'*_PSF_Library.pkl')[0]
        PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
        for i in range(len(PSF_list_clean)):
            psf = PSF_list_clean[i]
            data_process.PSF_list = [psf]
            pickle.dump(data_process , open('fit_material/'+'data_process_idx{2}_{0}_psf{1}.pkl'.format(filt, i, idx), 'wb'))
        
    for data_process in data_process_list_s:
        data_process.apertures = com_aper_s
        filt = data_process.filt
        PSF_lib_files = glob.glob('material/*'+filt[:-1]+'*_PSF_Library.pkl')[0]
        PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
        for i in range(len(PSF_list_clean)):
            psf = PSF_list_clean[i]
            data_process.PSF_list = [psf]
            pickle.dump(data_process , open('fit_material/'+'data_process_idx{2}_{0}_psf{1}.pkl'.format(filt, i, idx), 'wb'))
            