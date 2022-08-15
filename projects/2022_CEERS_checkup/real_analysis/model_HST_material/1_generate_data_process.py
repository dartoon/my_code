#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 23:13:44 2022

@author: Dartoon

Similar to JWST's 1_generate_data_process.py
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle
from galight.data_process import DataProcess
from galight.tools.plot_tools import plot_data_apertures_point
from galight.tools.cutout_tools import common_data_class_aperture
from galight.tools.measure_tools import measure_bkg
import warnings
warnings.filterwarnings("ignore")

shift_list = pickle.load(open('../model_JWST_stage3_Masafusa/material/jwst_shift_list.pkl','rb')) #For fov shift
ignore_id = [10, 21, 30, 31,  41, 46, 47, 52]
remove_id = [24, 55]

f = open("../model_JWST_stage3_Masafusa/material/target_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
result_folder = 'fit_result/'
lines = lines[1:]


HST_folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_HST_data/'
# for idx in range(len(lines)):
wht_max = {'F105W':11948740000, 'F125W':45840350000, 'F140W':858202000, 'F160W':45840350000.0}

for idx in range(1, 2):
    if idx in remove_id:
        continue
    line = lines[idx]
    HST_all_files= glob.glob(HST_folder+'/egs_all_wfc3_ir_*_030mas_v1.9_drz.fits')  #For NIRCam
    target_id, _RA, _Dec, spec_z, photo_z = line.split(' ')
    #load HST
    JWST_data_process_list = pickle.load(open('../model_JWST_stage3_Masafusa/material/data_process+apertures_{0}.pkl'.format(idx),'rb'))
    RA, Dec = JWST_data_process_list[0][0][0].cut_RA, JWST_data_process_list[0][0][0].cut_Dec, 
    com_aper_l = JWST_data_process_list[0][1]
    data_process_list = []
    # i = -1
    for i in range(len(HST_all_files)):
        HST_fitsFile = pyfits.open(HST_all_files[i])
        fov_image_HST = HST_fitsFile[0].data 
        header_HST = HST_fitsFile[0].header 
        if JWST_data_process_list[0][0] != []:
            exppix = 1
            radius = JWST_data_process_list[0][0][0].radius
        elif JWST_data_process_list[0][0] == []:
            exppix = 2
            radius = JWST_data_process_list[1][0][0].radius/1.5
        wht = pyfits.open(HST_all_files[i].replace('drz','wht'))[0].data
        exp = header_HST['EXPTIME']
        exp_map = exp * wht/ wht_max[header_HST['filter']]
        
        data_process = DataProcess(fov_image = fov_image_HST, target_pos = [RA, Dec], pos_type = 'wcs', header = header_HST,
                                  rm_bkglight = False, if_plot=False, zp = 27 , exptime= exp_map) #, fov_noise_map=fov_noise_map)
        #estimate local bkg and remove:
        data_process.generate_target_materials(radius=250, create_mask = False, npixels = 80,
                                                cut_kernel = None, skip = True)
        bkg_std = data_process.bkg_std
        bkglight = measure_bkg(data_process.target_stamp, if_plot=True) # measure local bkg light to remove
        
        data_process.target_pos = data_process.target_pos + np.array(shift_list[idx][0][-1]) * exppix # Match to JWST
        data_process.generate_target_materials(radius=radius,cut_kernel=None, create_mask = False, skip = True)
        ct = int((len(bkglight) - len(data_process.target_stamp ))/2)
        data_process.target_stamp = data_process.target_stamp - bkglight[ct:-ct, ct:-ct]
        
        data_process.target_stamp = np.flip(data_process.target_stamp)
        data_process.noise_map = np.flip(data_process.noise_map)
        if JWST_data_process_list[0][0] != []:
            data_process.target_mask = JWST_data_process_list[0][0][0].target_mask
        del data_process.fov_image
        del data_process.exptime
        del data_process.exp_map
        data_process.filt = header_HST['filter']
        if np.sum(data_process.target_stamp) != 0:
            data_process_list.append(data_process)
    filters =[data_process_list[i].filt for i in range(len(data_process_list))]
    if com_aper_l == []:
        com_aper_l = common_data_class_aperture(data_process_list, l_idx=0, return_idx=0)
    print("Common aperture:")
    for i in range(len(data_process_list)):
        print(idx, target_id, filters[i],':')
        plot_data_apertures_point(data_process_list[i].target_stamp * data_process_list[i].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                                  com_aper_l, figsize=(4,3))
    print("Above are for the", target_id, 'filts:', 
          filters, 'idx:', idx)
    hold = input('Hold ... OK?\n')
    pickle.dump([[data_process_list, com_aper_l]], open('material/'+'data_process+apertures_{0}.pkl'.format(idx), 'wb'))