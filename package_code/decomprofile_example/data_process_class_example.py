#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 16:01:17 2020

@author: Xuheng Ding

You can skip this step if the QSO stamp, noise level and the PSF is ready.
"""
#photutils in version 0.7.2
#astropy in version astropy-4.0.1


import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from matplotlib.colors import LogNorm

fitsFile = pyfits.open('../example_data/HST/QSO/1147_final_drz.fits')
fov_image = fitsFile[1].data # check the back grounp
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']

wht = fitsFile[2].data # The WHT map
import decomprofile.tools_data.astro_tools as astro_tools
exp =  astro_tools.read_fits_exp(fitsFile[0].header)  #Read the exposure time 
mean_wht = exp * (0.0642/0.135)**2
exp_map = exp * wht/mean_wht

from decomprofile.data_process import DataProcess

data_class = DataProcess(fov_image = fov_image, target_pos = [1135, 648], header = header, rm_bkglight = False, exptime = exp_map, if_plot=False)
data_class.generate_target_materials(target_radius=60, create_mask = True, nsigma=2.8, exp_sz= 1.2, npixels = 15, if_plot=True)
#%%
data_class.find_PSF(50)
# from decomprofile.tools_data.astro_tools import plt_fits
