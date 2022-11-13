#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 17:26:40 2022

@author: Dartoon
"""


import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle
import warnings
warnings.filterwarnings("ignore")

filt = 'F150W'
idx = 0
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# folder = '/Users/Dartoon/Downloads/z6JWSTNIRcam/NIRCam_{0}_{1}_stage2'.format(target_id[:5], filt)
folder = '/Users/Dartoon/Downloads/z6JWSTNIRcam/NIRCam_J2255_stage3_usehalf'
# RA, Dec = [[343.9071242140829, 2.8544273992380655],
#  [343.9087017408779, 2.873915874029475],
#  [343.91533911665135, 2.8554187194448932],
#  [343.897422888621, 2.8616087847280665],
#  [343.92142375403773, 2.8675199927646666],
#  [343.8991757740048, 2.8728876497240563],
#  [343.89762679815607, 2.8632665200635237],
#  [343.8925089271645, 2.873924069346915],
#  [343.9211321173894, 2.8462168153246497]]
filter_files= glob.glob(folder+'/*.fits'.format(target_id[:5], filt))  #For NIRCam
filter_files.sort()

#%%
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.tools.astro_tools import plt_fits
# filename = filter_files[1]
for filename in filter_files[:]:
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
    data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
                              rm_bkglight = True, exptime = exp_map, if_plot=False, zp = zp)#, fov_noise_map = fov_noise_map)
    # fitsFile[1].data = data_process.fov_image
    # fitsFile.writeto(folder+'/bkg_removed/'+filename.split('/')[-1][:-5] + '_rmbkg'+'.fits')
    
    data_process.fov_image = np.nan_to_num(data_process.fov_image)
    data_process.generate_target_materials(radius=30, create_mask = False, nsigma=1.5, 
                                            cut_kernel = None, if_select_obj=False,
                                          exp_sz= 1.2, npixels = 20, if_plot=False)
    
    plt_fits(data_process.target_stamp)