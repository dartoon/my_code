#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 10:42:55 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle
filters =  ['F150','F356']

folder = '../NIRCam_data/Dec7/'

filt = filters[0] #!!!

idx = 2

import sys
sys.path.insert(0, '../model_z6_data_id0/')
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

filter_files= glob.glob(folder+'/*{0}*{1}*i2d.fits'.format(target_id[:5], filt))  #For NIRCam
filter_files.sort()

folder = filter_files[0].split('SHELLQs_')[0]

#%%
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
# filename = filter_files[1]
for filename in filter_files:
    print("Select for", filename.split('/')[-1])
    # Grab the JWST provided ERR map:
    
    fitsFile = pyfits.open(filename)
    fov_image = fitsFile[1].data # check the back grounp
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        
    flux_mjsr = header['PHOTMJSR']
    pixscale = read_pixel_scale(header)
    zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
    
    # fov_noise_map = fitsFile[2].data
    
    wht = fitsFile[4].data # The WHT map
    exp = fitsFile[0].header['EFFEXPTM']
    gain_value = 2
    exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
        
    from galight.data_process import DataProcess
    data_process = DataProcess(fov_image = fov_image, target_pos = [600,1150], pos_type = 'pixel', header = header,
                              rm_bkglight = True, exptime = exp_map, if_plot=True, zp = zp)#, fov_noise_map = fov_noise_map)
    fitsFile[1].data = data_process.fov_image
    fitsFile.writeto(folder+'/bkg_removed/'+filename.split('/')[-1][:-5] + '_rmbkg'+'.fits')