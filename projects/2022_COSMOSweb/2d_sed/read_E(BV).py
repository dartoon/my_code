#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 14:19:19 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle, copy, matplotlib
from galight.tools.astro_tools import plt_fits_color
from photutils.aperture import RectangularAperture
from matplotlib.colors import LogNorm   

name = ['cid_473', 'cid_1210', 'cid_1245'][2]

image_list = pickle.load(open('colorimage_bin2_{0}.pkl'.format(name),'rb'))
images = []
# zp_list = []
color_list = [-1,-3,-4]
for i in color_list:  #['F115W', 'F150W','F277W', 'F444W']
    images.append(image_list[i])
from galight.tools.astro_tools import plt_fits_color, plt_fits
plt_fits_color(images, Q=7, stretch=0.3)

Ebv_image = pickle.load(open('E(BV)_{0}.pkl'.format(name),'rb'))
norm = None
# plt.imshow(Ebv_image, norm=norm, origin='lower' ) 
# plt.colorbar()
# plt.show()

f_center = len(Ebv_image)/2
deltaPix = 0.06  #Bin the 0.03 by factor of 2

if name == 'cid_473':
    w = 0.2 / deltaPix
    h = 1 / deltaPix
    theta = 20
if name == 'cid_1210':
    w = 0.2 / deltaPix
    h = 1 / deltaPix
    theta = -20
if name == 'cid_1245':
    w = 0.2 / deltaPix
    h = 1.5 / deltaPix
    theta = -20

aper = RectangularAperture((f_center, f_center), w, h, theta=theta/180*np.pi)
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat'))  # copy the default cmap
my_cmap.set_bad('black')
vmin = 1.e-3
# plt.imshow(Ebv_image, origin='lower', cmap=my_cmap, norm=LogNorm())  # , vmin=vmin, vmax=vmax)
plt.imshow(Ebv_image, norm=norm, origin='lower',vmin = vmin) 
plt.colorbar()
aper.plot(color='blue',
          lw=3.5)
plt.show()
from photutils.aperture import aperture_photometry
mask_ = aper.to_mask(method = 'center')
w,h = mask_.shape
mask = np.zeros_like(Ebv_image) 
mask[int(f_center - w/2): int(f_center + w/2) , int(f_center- h/2): int(f_center+h/2) ] +=  mask_
# plt.imshow(mask, origin='lower')
# plt.show()
mean = np.mean(Ebv_image[mask==1])
std = np.std(Ebv_image[mask==1])
print('{2} E(B-V) in aperture: {0:.3}$\pm${1:.3}'.format(mean,std, name))
# phot_table_host = to_mask(Ebv_image, aper)