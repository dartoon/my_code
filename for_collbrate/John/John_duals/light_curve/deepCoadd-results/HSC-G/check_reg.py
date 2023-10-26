#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 09:33:30 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt


import glob
files = glob.glob('8524/6,7/det*fits')

fitsFile = pyfits.open('/Users/Dartoon/Downloads/HSC_rerun/calexp-HSC-G-8524-6,7.fits')

plt.imshow(fitsFile[1].data)
plt.show()
# from galight.data_process import DataProcess



# #Data is available at: https://drive.google.com/file/d/1ZO9-HzV8K60ijYWK98jGoSoZHjIGW5Lc/view?usp=sharing
# fitsFile = pyfits.open('../example_files/HSC/QSO/000017.88+002612.6_HSC-I.fits')

# #Load the fov image data:
# fov_image = fitsFile[1].data # check the back grounp

# #Derive the header informaion, might be used to obtain the pixel scale and the exposure time.
# header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']

# #Derive the fov noise level map:
# err_data= fitsFile[3].data ** 0.5

# #Calculate the zeropoint for HSC filters:
# file_header0 = fitsFile[0].header
# FLUXMAG0 = file_header0['FLUXMAG0']
# zp =  2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.

# #Load the PSF data:
# PSF = pyfits.getdata('../example_files/HSC/QSO/000017.88+002612.6_HSC-I_psf.fits')

# #RA, DEC information of the QSO:
# QSO_RA, QSO_DEC = 0.07452999800443649, 0.4368380010128021
# data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [QSO_RA, QSO_DEC],
#                            pos_type = 'wcs', header = header,
#                           rm_bkglight = True, if_plot=False, zp = zp)

# #Generate the fitting materials
# data_process.generate_target_materials(radius=None, create_mask = True, nsigma=2.8,
#                                       exp_sz= 1.5, npixels = 15, if_plot=True)