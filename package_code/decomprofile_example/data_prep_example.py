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

#Load data and plot:
fitsFile = pyfits.open('../example_data/HST/QSO/1104_final_drz.fits')
img = fitsFile[1].data # check the back grounp
fig=plt.figure(figsize=(15,15))
ax=fig.add_subplot(1,1,1)
ax.imshow(img, norm=LogNorm(), origin='lower') 
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
plt.show()  
#%%
from astropy.visualization import SqrtStretch
from astropy.stats import SigmaClip
from photutils import Background2D, SExtractorBackground  
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import make_source_mask

norm = ImageNormalize(stretch=SqrtStretch())         
sigma_clip = SigmaClip(sigma=3., maxiters=10)
bkg_estimator = SExtractorBackground()
#Define the regions where contains the signal.
mask_0 = make_source_mask(img, nsigma=2, npixels=25, dilate_size=11) 
fig=plt.figure(figsize=(15,15))
ax=fig.add_subplot(1,1,1)
ax.imshow(mask_0, origin='lower') 
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
plt.show() 

#%%
mask_1 = (np.isnan(img))
mask = mask_0 + mask_1
#estimate the 2D background light:
bkg = Background2D(img, (50, 50), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                   mask=mask)

back = bkg.background* ~mask_1  #The 2-D back ground light estimated from the whole field.
fig=plt.figure(figsize=(15,15))
ax=fig.add_subplot(1,1,1)
ax.imshow(back, origin='lower', cmap='Greys_r')
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
plt.show()

img = img - back              
# pyfits.PrimaryHDU(img).writeto('image_subbkg.fits',overwrite=True)

#%%
from decomprofile.tools_data.image_tailor import cut_center_bright
QSO_loc = [1138, 648]  # The postion of the QSO in the frame
QSO, cut_center = cut_center_bright(image=img, center= QSO_loc,  kernel = 'center_gaussian', radius=60, return_center=True, if_plot=True)
[] TODO: Test bright center and Gussian Center
