#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 09:44:37 2018

@author: Dartoon

Test Background remove.
"""
from photutils.datasets import make_100gaussians_image
data = make_100gaussians_image()
import matplotlib.pyplot as plt

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
norm = ImageNormalize(stretch=SqrtStretch())
plt.imshow(data, norm=norm, origin='lower', cmap='Greys_r')
plt.show()

import numpy as np
from astropy.stats import biweight_location
print(np.median(data))
print(biweight_location(data))
from astropy.stats import mad_std
print(mad_std(data))

from astropy.stats import sigma_clipped_stats
mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)
print((mean, median, std))

from photutils import make_source_mask
mask = make_source_mask(data, snr=2, npixels=5, dilate_size=11)
mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
print((mean, median, std))  

ny, nx = data.shape
y, x = np.mgrid[:ny, :nx]
gradient =  x * y / 5000.
data2 = data + gradient
plt.imshow(data2, norm=norm, origin='lower', cmap='Greys_r') 
plt.show()

from astropy.stats import SigmaClip
from photutils import Background2D, SExtractorBackground
sigma_clip = SigmaClip(sigma=3., iters=10)
bkg_estimator = SExtractorBackground()
bkg = Background2D(data2, (50, 50), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
plt.imshow(data2, norm=norm, origin='lower', cmap='Greys_r') 
bkg.plot_meshes(outlines=True, color='#1f77b4')
plt.show()

print(bkg.background_median)
print(bkg.background_rms_median)
plt.imshow(bkg.background, origin='lower', cmap='Greys_r')
plt.show()

#from scipy.ndimage import rotate
#data3 = rotate(data2, -45.)
#norm = ImageNormalize(stretch=SqrtStretch())    
#plt.imshow(data3, origin='lower', cmap='Greys_r', norm=norm) 
#plt.show()   
#mask = (data3 == 0)
#bkg3 = Background2D(data3, (25, 25), filter_size=(3, 3), mask=mask)
#
#back3 = bkg3.background * ~mask
#norm = ImageNormalize(stretch=SqrtStretch())    
#plt.imshow(back3, origin='lower', cmap='Greys_r', norm=norm)   
#plt.show()
#plt.imshow(data3, origin='lower', cmap='Greys_r', norm=norm)
#bkg3.plot_meshes(outlines=True, color='#1f77b4')
#plt.show()