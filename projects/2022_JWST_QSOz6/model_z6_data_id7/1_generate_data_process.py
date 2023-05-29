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


data_type = 'all' 
filt = 'F150W'
file_NO = 0

fov = 'small'

idx = 7
# folder = '/Users/Dartoon/Downloads/z6JWSTNIRcam/NIRCam_J2255_stage3_{0}/bkg_removed'.format(data_type)
folder = '../NIRCam_data/Apr27/bkg_removed/' 
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# jwst_all_filenames = glob.glob(folder+'/*{0}*{1}*.fits'.format(target_id[:5], filts[0]))
jwst_all_filenames = glob.glob(folder+'*{0}*{1}*_rmbkg.fits'.format(target_id[:5], filt))  #For NIRCam
jwst_all_filenames.sort()
file = jwst_all_filenames[file_NO]
if data_type == 'all':
    if fov != 'large':
        run_folder = 'stage3_{0}/'.format(data_type)
    else:
        run_folder = 'stage3_{0}_largecut/'.format(data_type)
    
elif data_type == 'half':
    if file_NO == 0:
        run_folder = 'stage3_first_half/'
    if file_NO == 1:
        run_folder = 'stage3_second_half/'

result_folder = run_folder + 'fit_result/'

if filt == 'F356W':
    if fov == 'large':
        radius = 80
    else:
        radius = 35
elif filt == 'F150W':
    if fov == 'large':
        radius = 40
    else:
        radius = 35


#%%
cut_kernel = None #After pos correct then, do the nearest_obj_center
# filts = ['F356W', 'F150W']
for filt in [filt]:
    _fitsFile = pyfits.open(file)
    fov_image = _fitsFile[1].data # check the back grounp
    header = _fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    flux_mjsr = header['PHOTMJSR']
    pixscale = read_pixel_scale(header)
    zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
    wht = _fitsFile[4].data # The WHT map
    exp = _fitsFile[0].header['EFFEXPTM']
    print("Exp time:", exp)
    if _fitsFile[0].header['CHANNEL'] == 'LONG':
        gain_value = 2
        expsize = 1
    else:
        gain_value = 1.8
        expsize = 1
    exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
    print("Processing data...")
    data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], #The final cut center is dete
                                   pos_type = 'wcs', header = header,rm_bkglight = False, 
                                   if_plot=False, zp = zp, exptime= exp_map, 
                                   fov_noise_map = None)
    
    data_process.generate_target_materials(radius=radius * expsize, create_mask = False, nsigma=1.5, 
                                            cut_kernel = 'center_bright', if_select_obj=False,
                                          exp_sz= 1.2, npixels = 60 * expsize, if_plot=False)
    
    # # data_process.apertures = []
    # data_process.apertures[0].theta = 0
    # data_process.apertures[0].b = data_process.apertures[0].b/2
    # if _fitsFile[0].header['CHANNEL'] == 'LONG':
    #     data_process.apertures[0].positions[0] = data_process.apertures[0].positions[0] - 1.5
    # if _fitsFile[0].header['CHANNEL'] != 'LONG':
    #     data_process.apertures[0].positions[0] = data_process.apertures[0].positions[0] - 1.5/2
    del data_process.fov_image
    del data_process.exptime
    
    #%%
    
    # ap = data_process.apertures[2]
    # data_process.apertures[2] = data_process.apertures[3]
    # data_process.apertures[3] = ap
    print(target_id, filt, 'apertures', len(data_process.apertures) )
    data_process.filt = filt
    data_process.file = file
    data_process.plot_aperture()
        
    hold = input('Hold ... OK?\n')
    # if hold == '!':
    #     break
    pickle.dump(data_process, open(run_folder+'material/'+'data_process_idx{0}_{1}.pkl'.format(idx, filt), 'wb'))

# import shutil
# shutil.copyfile('1_generate_data_process.py', 'backup_files/1_generate_data_process_idx{0}.py'.format(idx))