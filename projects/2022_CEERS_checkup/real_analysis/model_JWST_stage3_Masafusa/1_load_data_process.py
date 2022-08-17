#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 15:01:12 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from galight.tools.plot_tools import plot_data_apertures_point
import pickle
import glob

f = open("material/target_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
result_folder = 'fit_result/'
lines = lines[1:]
remove_id = [24, 55]
cid = 0

# files = glob.glob('material/data_process+apertures_*.pkl')
# files.sort()
for idx in range(cid,cid+10):
    if idx in remove_id:
        continue
    line = lines[idx]
    target_id, RA, Dec, spec_z, photo_z = line.split(' ')
    file = 'material/data_process+apertures_{0}.pkl'.format(idx)
    data_process_aperture_list = pickle.load(open(file,'rb')) #For fov shift
    [[data_process_list_l, com_aper_l], [data_process_list_s, com_aper_s] ] = data_process_aperture_list
    
    print("Common aperture:")
    for i in range(len(data_process_list_l)):
        print(idx, target_id, data_process_list_l[i].filt,':')
        plot_data_apertures_point(data_process_list_l[i].target_stamp * data_process_list_l[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                                  com_aper_l, figsize=(4,3))
    for i in range(len(data_process_list_s)):
        print(idx, target_id, data_process_list_s[i].filt,':')
        plot_data_apertures_point(data_process_list_s[i].target_stamp * data_process_list_s[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                                  com_aper_s, figsize=(4,3))