#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 15:49:05 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from galight.tools.plot_tools import plot_data_apertures_point
import pickle
import glob
from astropy.wcs import WCS

f = open("material/target_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
result_folder = 'fit_result/'
lines = lines[1:]
remove_id = [24, 55]
cid = 0
import warnings
warnings.filterwarnings("ignore")

# files = glob.glob('material/data_process+apertures_*.pkl')
# files.sort()
for idx in [0, 1, 2, 28, 35, 51]:  #z_spec > 1.6
# for idx in range(9, len(lines)):
# for idx in [1]:
    if idx in remove_id:
        continue
    line = lines[idx]
    target_id, RA, Dec, spec_z, photo_z = line.split(' ')
    HST_file = '../model_HST_material/material/data_process+apertures_{0}_IR.pkl'.format(idx)
    HST_data_process_aperture_list = pickle.load(open(HST_file,'rb')) #For fov shift
    [[data_process_list_HST, com_aper_HST]] = HST_data_process_aperture_list
    HST_ACS_file = '../model_HST_material/material/data_process+apertures_{0}_ACS.pkl'.format(idx)
    HST_data_process_aperture_list_ACS = pickle.load(open(HST_ACS_file,'rb')) #For fov shift
    [[data_process_list_HST_ACS, com_aper_HST_ACS]] = HST_data_process_aperture_list_ACS
      
    file = 'material/data_process+apertures_{0}.pkl'.format(idx)
    data_process_aperture_list = pickle.load(open(file,'rb')) #For fov shift
    [[data_process_list_l, com_aper_l], [data_process_list_s, com_aper_s] ] = data_process_aperture_list
    
    print("Common aperture:", RA, Dec)

    for i in range(len(data_process_list_HST_ACS)):
        print(idx, target_id, data_process_list_HST_ACS[i].filt,':')
        plot_data_apertures_point(data_process_list_HST_ACS[i].target_stamp * data_process_list_HST_ACS[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                                  com_aper_HST_ACS, figsize=(4,3))
        
    for i in range(len(data_process_list_HST)):
        print(idx, target_id, data_process_list_HST[i].filt,':')
        plot_data_apertures_point(data_process_list_HST[i].target_stamp * data_process_list_HST[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                                  com_aper_HST, figsize=(4,3))
    
    for i in range(len(data_process_list_l)):
        print(idx, target_id, data_process_list_l[i].filt,':')
        plot_data_apertures_point(data_process_list_l[i].target_stamp * data_process_list_l[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                                  com_aper_l, figsize=(4,3))
    for i in range(len(data_process_list_s)):
        print(idx, target_id, data_process_list_s[i].filt,':')
        plot_data_apertures_point(data_process_list_s[i].target_stamp * data_process_list_s[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                                  com_aper_s, figsize=(4,3))
    
    for i in range(len(data_process_list_HST_ACS)):
        _data_process = data_process_list_HST_ACS[i]
        header = _data_process.header
        wcs = WCS(header)
        cut_RA, cut_Dec = wcs.all_pix2world([_data_process.target_pos], 1)[0] #re-define RA, Dec
        shift = (cut_RA - float(RA))*3600/_data_process.deltaPix , (cut_Dec-float(Dec))*3600/_data_process.deltaPix
        print('HST shift', data_process_list_HST_ACS[i].filt, shift)
        
    for i in range(len(data_process_list_HST)):
        _data_process = data_process_list_HST[i]
        header = _data_process.header
        wcs = WCS(header)
        cut_RA, cut_Dec = wcs.all_pix2world([_data_process.target_pos], 1)[0] #re-define RA, Dec
        shift = (cut_RA - float(RA))*3600/_data_process.deltaPix , (cut_Dec-float(Dec))*3600/_data_process.deltaPix
        print('HST shift', data_process_list_HST[i].filt, shift)
    
    for i in range(len(data_process_list_l)):
        _data_process = data_process_list_l[i]
        header = _data_process.header
        wcs = WCS(header)
        cut_RA, cut_Dec = wcs.all_pix2world([_data_process.target_pos], 1)[0] #re-define RA, Dec
        shift = (cut_RA - float(RA))*3600/_data_process.deltaPix , (cut_Dec-float(Dec))*3600/_data_process.deltaPix
        print('LW shift', data_process_list_l[i].filt, shift)
        
    for i in range(len(data_process_list_s)):
        _data_process = data_process_list_s[i]
        header = _data_process.header
        wcs = WCS(header)
        cut_RA, cut_Dec = wcs.all_pix2world([_data_process.target_pos], 1)[0] #re-define RA, Dec
        shift = (cut_RA - float(RA))*3600/_data_process.deltaPix , (cut_Dec-float(Dec))*3600/_data_process.deltaPix
        print('SW shift',data_process_list_s[i].filt , shift)
        
    hold = input(str(idx)+', idx OK?')