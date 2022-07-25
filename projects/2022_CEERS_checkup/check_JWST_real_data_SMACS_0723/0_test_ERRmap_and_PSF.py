#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 13:25:08 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale

#%%
fitsFile = pyfits.open('/Users/Dartoon/Downloads/JWST_ID2736/F356W/jw02736-o001_t001_nircam_clear-f356w/jw02736-o001_t001_nircam_clear-f356w_i2d.fits')
fov_image = fitsFile[1].data # check the back grounp
fov_noise_image = fitsFile[2].data # check the back grounp
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']

flux_mjsr = header['PHOTMJSR']
pixscale = read_pixel_scale(header)
zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux

#%% Grab the JWST provided ERR map:
data_process = DataProcess(fov_image = fov_image, target_pos = [700., 500.], pos_type = 'pixel', header = header, fov_noise_map = fov_noise_image,
                          rm_bkglight = True, if_plot=False, zp = zp)
data_process.generate_target_materials(radius=65, create_mask = False, nsigma=2.8, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 15, if_plot=True)
noise_map_ERR = data_process.noise_map

#%% Calculated my Error map:
wht = fitsFile[4].data # The WHT map
exp =  header['XPOSURE']  #Read the exposure time 
gain_value = 2
exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
data_process = DataProcess(fov_image = fov_image, target_pos = [700., 500.], pos_type = 'pixel', header = header,
                          rm_bkglight = True, if_plot=False, zp = zp, exptime= exp_map )
data_process.generate_target_materials(radius=65, create_mask = False, nsigma=2.8, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 15, if_plot=True)
noise_map_cal = data_process.noise_map
# pyfits.PrimaryHDU(noise_map_ERR).writeto('local_error_JWST.fits',overwrite=True)
# pyfits.PrimaryHDU(noise_map_cal).writeto('local_error_galight.fits',overwrite=True)

#%%PSF works.
data_process.find_PSF(radius = 40, user_option = True, if_filter=True, FWHM_filer=2.35, neighborhood_size=50,
                      nearyby_obj_filter=True, FWHM_sort=True)
# data_process.find_PSF(radius = 50, PSF_pos_list = [[ 350., 1055.], [2078., 1910.]], user_option = False)
data_process.plot_overview(label = 'Example', target_label = None)

# from galight.tools.cutout_tools import psf_clean
# psf_clean(data_process.PSF_list[3])