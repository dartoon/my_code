#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 15:49:55 2018

@author: Dartoon

find the local maximum as point position
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

data = pyfits.getdata('test_image.fits')
neighborhood_size = 4
threshold = 5

# filters.maximum_filter is to make the pixel value as the value as the local max
data_max = filters.maximum_filter(data, neighborhood_size) 
maxima = (data == data_max)
data_min = filters.minimum_filter(data, neighborhood_size)
diff = ((data_max - data_min) > threshold)
maxima[diff == 0] = 0
#plt.imshow(data_max, origin='low')
#plt.colorbar()
#plt.show()
#
#plt.imshow(data_min, origin='low')
#plt.colorbar()
#plt.show()
      
labeled, num_objects = ndimage.label(maxima)
slices = ndimage.find_objects(labeled)

x, y = [], []
for dy,dx in slices:
    x_center = (dx.start + dx.stop - 1)/2
    x.append(x_center)
    y_center = (dy.start + dy.stop - 1)/2    
    y.append(y_center)

plt.imshow(data, origin='low')
plt.colorbar()
plt.plot(x,y, 'ro')
plt.show()