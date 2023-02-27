#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 14:47:09 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

#Plot 1D:
# SDSS2304-0038
folder = 'MAST_2023-02-21T0308/HST/iexo11020/'
filename = 'iexo11020_drz.fits'
fitsFile = pyfits.open(folder+filename)
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
img = fitsFile[1].data #
img_sub = img[561:568,610:750]
from galight.tools.astro_tools import plt_fits
plt_fits(img_sub,figsize=(11, 8))
spec_1 = img_sub[4:,:]
spec_2 = img_sub[0:3,:]
z = 2.775
MgII_pos = 2798 * (1+z)
# CIII_pos = 1909 * (1+z)
shift = MgII_pos - 2500
plt_fits(spec_1, figsize=(11, 8))
plt.figure(figsize=(11, 8))
plt.plot( (np.arange(0, len(spec_1.T))*24.5 + shift),np.sum(spec_1,axis=0))
plt.xlabel('$\AA$')
plt.show()
plt_fits(spec_1, figsize=(11, 8))
plt.figure(figsize=(11, 8))
plt.plot( (np.arange(0, len(spec_1.T))*24.5 + shift)/(1+z),np.sum(spec_1,axis=0))
plt.xlabel('rest-frame  $\AA$')
plt.show()
plt_fits(spec_2,figsize=(11, 8))
plt.figure(figsize=(11, 8))
plt.plot( (np.arange(0, len(spec_2.T))*24.5 + shift),np.sum(spec_2,axis=0))
plt.xlabel('$\AA$')
plt.show()
plt_fits(spec_2,figsize=(11, 8))
plt.figure(figsize=(11, 8))
plt.plot( (np.arange(0, len(spec_2.T))*24.5 + shift)/(1+z),np.sum(spec_2,axis=0))
plt.xlabel('rest-frame $\AA$')
plt.show()

#%%
#SDSS1246-0017
folder = 'MAST_2023-02-23T2316/HST/iexo07020/'
filename = 'iexo07020_drz.fits'
fitsFile = pyfits.open(folder+filename)
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
img = fitsFile[1].data #
img_sub = img[561:570,580:750]
from galight.tools.astro_tools import plt_fits

plt_fits(img_sub, figsize=(11,8))
z = 2.559
MgII_pos = 2798 * (1+z)
# CIII_pos = 1909 * (1+z)
shift = MgII_pos - 2600
w = 170
h = 3.5
theta = 0.6
from photutils.aperture import RectangularAperture
aper = RectangularAperture((85,5), w, h, theta=theta/180*np.pi)
mask_ = aper.to_mask(method = 'exact')
mask_1 = mask_.to_image(img_sub.shape)
plt.figure(figsize=(11, 8))
spec_1 = img_sub*mask_1
plt.imshow(spec_1, norm=None, origin='lower',vmin = 0.01, cmap = 'gist_heat') 
aper.plot(color='blue',lw=2.5)
plt.show()
plt.figure(figsize=(11, 8))
plt.plot( (np.arange(2, len(spec_1.T))*24.5 + shift),np.sum(spec_1,axis=0)[2:])
plt.xlabel('$\AA$')
plt.show()
plt.figure(figsize=(11, 8))
plt.plot( (np.arange(2, len(spec_1.T))*24.5 + shift)/(1+z),np.sum(spec_1,axis=0)[2:])
plt.xlabel('rest-frame $\AA$')
plt.show()
aper = RectangularAperture((85,2), w, h, theta=theta/180*np.pi)
mask_ = aper.to_mask(method = 'exact')
mask_2 = mask_.to_image(img_sub.shape)
plt.figure(figsize=(11, 8))
spec_2 = img_sub*mask_2
plt.imshow(spec_2, norm=None, origin='lower',vmin = 0.01, cmap = 'gist_heat') 
aper.plot(color='blue',lw=2.5)
plt.show()
plt.figure(figsize=(11, 8))
plt.plot( (np.arange(5, len(spec_2.T))*24.5 + shift),np.sum(spec_2,axis=0)[5:])
plt.xlabel('$\AA$')
plt.show()
plt.figure(figsize=(11, 8))
plt.plot( (np.arange(5, len(spec_2.T))*24.5 + shift)/(1+z),np.sum(spec_2,axis=0)[5:])
plt.xlabel('rest-frame $\AA$')
plt.show()
