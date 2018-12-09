#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 18:24:48 2018

@author: Dartoon

test load data
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

fitsFile = pyfits.open('../IMAGES/2-cutout-HSC-G-9558-pdr1_wide.fits')
sci_data= fitsFile[1].data.T
                  
plt.imshow(sci_data.T, origin='low', norm= LogNorm())
plt.colorbar()
plt.show()

print sci_data[331, 141]

