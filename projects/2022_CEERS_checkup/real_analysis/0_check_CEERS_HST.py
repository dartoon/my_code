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
# =============================================================================
# # JWST example
# # folder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_CEERS/CEERS_JWST_data'
# # all_files= glob.glob(folder+'/*clear*/*_i2d.fits')  #For NIRCam
# # fitsFile = pyfits.open(all_files[0])
# # fov_image = fitsFile[1].data # check the back grounp
# # header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
# =============================================================================
folder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_CEERS/CEERS_HST_data'
all_files= glob.glob(folder+'/egs_all_wfc3_ir_f160w_030mas_v1.9_drz.fits')  #For NIRCam
fitsFile = pyfits.open(all_files[0])
fov_image = fitsFile[0].data # check the back grounp
header = fitsFile[0].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
from galight.tools.astro_tools import plt_fits
plt_fits(fov_image[14000:14000+8000,19000:19000+8000])
# from galight.data_process import DataProcess

# data_process = DataProcess(fov_image = fov_image, target_pos = [5000, 5000], pos_type = 'pixel', header = header,
#                           rm_bkglight = False, if_plot=False, zp = 27 )
# # data_process.generate_target_materials(radius=30, create_mask = False, nsigma=2.8, if_select_obj=False,
# #                                       exp_sz= 1.2, npixels = 15, if_plot=True)
# # 
# # #%%PSF works.
# data_process.find_PSF(radius = 40, user_option = True, if_filter=False, neighborhood_size=50,
#                       nearyby_obj_filter=True, FWHM_sort=True, select_all=False)
# # data_process.plot_overview(label = 'Example', target_label = None)