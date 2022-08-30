#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 00:26:46 2022

@author: Dartoon

Load the CEERS HST image for quick check
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# # JWST example
# # folder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_CEERS/CEERS_JWST_data'
# # all_files= glob.glob(folder+'/*clear*/*_i2d.fits')  #For NIRCam
# # fitsFile = pyfits.open(all_files[0])
# # fov_image = fitsFile[1].data # check the back grounp
# # header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
# =============================================================================
folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_HST_data/'
all_files= glob.glob(folder+'/egs_all_acs_wfc_f606w_030mas_v1.9_drz.fits')  #For NIRCam
fitsFile = pyfits.open(all_files[0])
fov_image = fitsFile[0].data # check the back grounp
header = fitsFile[0].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
from galight.tools.astro_tools import plt_fits
# plt_fits(fov_image[14000:14000+8000,19000:19000+8000])

wht = pyfits.open(all_files[0].replace('drz','wht'))[0].data

# HST_all_files= glob.glob(folder+'/egs_all_acs_wfc*_030mas_v1.9_drz.fits')  #For NIRCam

#%%
# target_id = 'aegis_533'
# RA, Dec = 214.87553, 52.866464

# from galight.data_process import DataProcess
# data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
#                           rm_bkglight = False, if_plot=False, zp = 27 )
# data_process.generate_target_materials(radius=60, create_mask = False, nsigma=2.8, if_select_obj=False,
#                                       exp_sz= 1.2, npixels = 15, if_plot=False, cut_kernel = None)
# data_process.plot_aperture()

