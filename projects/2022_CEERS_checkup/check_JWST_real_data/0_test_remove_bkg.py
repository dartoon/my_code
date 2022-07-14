#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:35:56 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

# /Users/Dartoon/Downloads/JWST_ID2736/F356W/jw02736-o001_t001_nircam_clear-f356w/jw02736-o001_t001_nircam_clear-f356w_i2d.fits


#Test remove the bkg of final drizzled image:
#%%
fitsFile = pyfits.open('/Users/Dartoon/Downloads/JWST_ID2736/F356W/jw02736-o001_t001_nircam_clear-f356w/jw02736-o001_t001_nircam_clear-f356w_i2d.fits')
# fitsFile = pyfits.open('/Users/Dartoon/Downloads/JWST_ID2736/F150W/jw02736-o001_t001_nircam_clear-f150w/jw02736-o001_t001_nircam_clear-f150w_i2d.fits')
# fitsFile = pyfits.open('/Users/Dartoon/Downloads/JWST_ID2736/F356W/jw02736001001_02103_00001_nrcalong/jw02736001001_02103_00001_nrcalong_i2d.fits')

fov_image = fitsFile[1].data # check the back grounp
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']

#
wht = fitsFile[2].data # The WHT map
import galight.tools.astro_tools as astro_tools
exp =   header['XPOSURE'] #Read the exposure time 
mean_wht = exp #* (0.0642/0.135)**2
exp_map = exp * wht/mean_wht

#Start to use galight
from galight.data_process import DataProcess
data_process = DataProcess(fov_image = fov_image, target_pos = [50., 50.], pos_type = 'pixel', header = header,
                          rm_bkglight = True, exptime = exp_map, if_plot=True, zp = 27.0)
# pyfits.PrimaryHDU(data_process.fov_image,header=header).writeto('bkg_remove.fits',overwrite=True)

#It works really well!