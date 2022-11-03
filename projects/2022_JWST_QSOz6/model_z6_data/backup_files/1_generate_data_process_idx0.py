#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 16:13:10 2022

@author: Dartoon

After 1_align_jwst_astromerty.py. 

Propose: Generate the data_process for JWST and HST together. Apertures will also be generated. 

Template: 1_test_model_allbands.py 
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import sys
import pickle
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.tools.measure_tools import measure_bkg
from astropy.wcs import WCS
from galight.tools.cutout_tools import common_data_class_aperture
from galight.tools.plot_tools import plot_data_apertures_point
from galight.tools.cutout_tools import cutout
import warnings
warnings.filterwarnings("ignore")

folder = '/Users/Dartoon/Downloads/z6JWSTNIRcam/bkg_removed'

idx = 0
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

filts = ['F150W']
jwst_all_filenames = glob.glob(folder+'/*{0}*{1}*.fits'.format(target_id[:5], filts[0]))
file = jwst_all_filenames[0]
result_folder = 'fit_result/'
#%%
cut_kernel = None #After pos correct then, do the nearest_obj_center
# filts = ['F356W', 'F150W']
for filt in filts:
    _fitsFile = pyfits.open(file)
    fov_image = _fitsFile[1].data # check the back grounp
    header = _fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    flux_mjsr = header['PHOTMJSR']
    pixscale = read_pixel_scale(header)
    zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
    wht = _fitsFile[4].data # The WHT map
    exp = _fitsFile[0].header['EFFEXPTM']
    if _fitsFile[0].header['CHANNEL'] == 'LONG':
        gain_value = 2
        expsize = 1
        exppix = 1
    else:
        gain_value = 1.8
        expsize = 1.4
        exppix = 2
    exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
    print("Processing data...")
    data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], #The final cut center is dete
                                   pos_type = 'wcs', header = header,rm_bkglight = False, 
                                   if_plot=False, zp = zp, exptime= exp_map, 
                                   fov_noise_map = None)
    
    data_process.generate_target_materials(radius=30 * expsize, create_mask = False, nsigma=1.5, 
                                            cut_kernel = None, if_select_obj=False,
                                          exp_sz= 1.2, npixels = 40 * expsize, if_plot=False)
    
    # data_process.apertures = []
    data_process.apertures[0].b = data_process.apertures[0].b/2
    data_process.apertures[0].positions[0] = data_process.apertures[0].positions[0] - 3
    
    del data_process.fov_image
    del data_process.exptime
    print(target_id, filt, 'apertures', len(data_process.apertures) )
    data_process.filt = filt
    data_process.file = file
    data_process.plot_aperture()
        
    hold = input('Hold ... OK?\n')
    if hold == '!':
        break
    pickle.dump(data_process, open('material/'+'data_process_idx{0}_{1}.pkl'.format(idx, filt), 'wb'))

import shutil
shutil.copyfile('1_generate_data_process.py', 'backup_files/1_generate_data_process_idx{0}.py'.format(idx))