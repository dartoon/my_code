#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 14:49:32 2024

@author: Dartoon
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

# =============================================================================
# Notes: idx 2 F150W has some WCS offset
# =============================================================================


import sys
sys.path.insert(0, '../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info

filt = 'F150W' #!!!

idx = 9 #!!!
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
    
folder = '../../2022_JWST_QSOz6/NIRCam_data/*/bkg_removed/'   #!!!
jwst_all_filenames = glob.glob(folder+'*{0}*{1}*{2}*_rmbkg.fits'.format(target_id[:5],target_id[-4:], filt))  #For NIRCam
jwst_all_filenames.sort()

file = jwst_all_filenames[0]
_fitsFile = pyfits.open(file)
# fov_image = _fitsFile[1].data # check the back grounp
header = _fitsFile[1].header
# print(read_pixel_scale(header))


# jwst_all_filenames = glob.glob(folder+'/*{0}*{1}*.fits'.format(target_id[:5], filts[0]))
jwst_all_filenames = glob.glob(folder+'*{0}*{1}*{2}*_rmbkg.fits'.format(target_id[:5],target_id[-4:], filt))  #For NIRCam
jwst_all_filenames.sort()
if len(jwst_all_filenames) > 1 :
    warnings.warn('Find more than one JWST file')
    

if filt == 'F356W':
    radius = 40
elif filt == 'F150W':
    radius = 35

# filts = ['F356W', 'F150W']
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
else:
    gain_value = 1.8
exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
print("Processing data...")
data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], #The final cut center is dete
                                pos_type = 'wcs', header = header,rm_bkglight = False, 
                                if_plot=False, zp = zp, exptime= exp_map, 
                                fov_noise_map = None)

data_process.generate_target_materials(radius=radius, create_mask = False, nsigma=1.5, 
                                        cut_kernel = 'center_bright', if_select_obj=False,
                                      exp_sz= 1.2, npixels = 50 , if_plot=False)

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

save_folder = '../'
hold = input('Hold ... OK?\n')
pickle.dump(data_process, open(save_folder+'material/{0}/'.format(filt)+'data_process_idx{0}.pkl'.format(idx), 'wb'))

# import shutil
# shutil.copyfile('1_generate_data_process.py', 'backup_files/1_generate_data_process_idx{0}.py'.format(idx))