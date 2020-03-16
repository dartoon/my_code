#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 15:18:14 2017

@author: dxh
"""
import numpy as np
import matplotlib.pyplot as plt

import astropy.io.fits as pyfits
from astrofunc.LightProfiles.shapelets import ShapeletSet
shapeletSet = ShapeletSet()
import astrofunc.util as util
from astrofunc.util import Util_class
util_class = Util_class()
#==============================================================================
# #CUT M31
#==============================================================================
from source_info import  source_list #[0]: file name [1]: total size [2]: galfit R_e [3]:R_e/totalzise
source_name=source_list[0][1]
name=source_name
print "source name:", name
hdu = pyfits.open('fix/{0}_fix.fits'.format(name))
data = hdu[0].data
hdu.close()
#data -= data.min()
plt.matshow(np.log10(data),origin='lower')
plt.colorbar()
plt.show()
data /= data.sum()/1000.

if len(data)>300:
    factor=2
    if len(data)%factor !=0:
        res=len(data)%factor
        data_g=data[:-res,:-res]
        data=rebin.block(data_g,(len(data_g)/factor,len(data_g)/factor),factor=factor)
    else:
        data=rebin.block(data,(len(data)/factor,len(data)/factor),factor=factor)

x, y = util.make_grid(numPix=len(data), deltapix=1)  # make a coordinate grid

image_1d = util.image2array(data)  # map 2d image in 1d data array

n_max = 100  # choice of number of shapelet basis functions, 150 is a high resolution number, but takes long
beta = len(data)/20.  # shapelet scale parameter
# return the shapelet coefficients
param_list = shapeletSet.decomposition(image_1d, x, y, n_max, beta, 1., center_x=0, center_y=0) 
# reconstruct M31 with the shapelet coefficients
image_reconstructed = shapeletSet.function(x, y, param_list, n_max, beta, center_x=0, center_y=0)
image_reconstructed_2d = util.array2image(image_reconstructed)  # map 1d data vector in 2d image
plt.matshow(np.log10(image_reconstructed_2d),origin='lower')#,vmin=-5, vmax=1)
plt.colorbar()
plt.show()
pyfits.PrimaryHDU(image_reconstructed_2d).writeto('shapelet/{0}_shapelet.fits'.format(name),overwrite=True)

import pickle
output= open('shapelet/{0}.pkl'.format(name),'wb')
passvar={'n_max': n_max, 'beta': beta, 'param_list': param_list, 'numPix_large': len(data)}
pickle.dump(passvar,output)