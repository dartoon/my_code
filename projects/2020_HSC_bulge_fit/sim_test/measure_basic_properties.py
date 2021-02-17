#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 22:36:12 2021

@author: Dartoon
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from decomprofile.tools.measure_tools import measure_FWHM
import glob

files = glob.glob("../SDSS_0.2-0.3/*_psf.fits")

# FWHMs = []
# for file in files:
#     image = pyfits.getdata(file)
#     # print(image.shape)
#     FWHM = np.mean( measure_FWHM(image) ) * 0.165
#     FWHMs.append(FWHM)

# #%%    
# plt.hist(FWHMs)
# plt.xlabel('arcsec')
# plt.show()

#%%
# from decomprofile.tools.measure_tools import esti_bgkstd
# from decomprofile.tools.measure_tools import measure_bkg
# stds = []
# # for i in range(len(files)):
# # for file in files:
# for i in [427]:     
#     file = files[i]
#     file = file.split('_psf')[0]+file.split('_psf')[1]
#     image = pyfits.getdata(file)
#     bkglight = measure_bkg(image)
#     image = image-bkglight
#     std = esti_bgkstd(image[50:-50,50:-50],  nsigma=2, exp_sz= 1.5, npixels = 15, if_plot=False)
#     stds.append(std)    
# plt.hist(stds)
# plt.xlabel('std')
# plt.show()

#%%
from decomprofile.tools.astro_tools import read_fits_exp

exps = []
for file in files:
    file = file.split('_psf')[0]+file.split('_psf')[1]
    fitsFile = pyfits.open(file)
    header = fitsFile[1].header 
    # exp = read_fits_exp(header)
    # # print(image.shape)
    # exps.append(exp)
