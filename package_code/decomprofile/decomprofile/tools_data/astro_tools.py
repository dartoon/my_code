#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 14:32:21 2020

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

from astropy.wcs import WCS
from matplotlib.colors import LogNorm

def read_pixel_scale(fitsFile,frame=0):
    wcs = WCS(fitsFile[frame].header)
    diff_RA_DEC = wcs.all_pix2world([0,0],[0,1],1)
    diff_scale = np.sqrt((diff_RA_DEC[0][1]-diff_RA_DEC[0][0])**2 + (diff_RA_DEC[1][1]-diff_RA_DEC[1][0])**2)
    pix_scale = diff_scale * 3600
    return pix_scale

def read_fits_exp(fitsFile,frame=0):
    file_header0 = fitsFile[frame].header
    return file_header0['EXPTIME']

def plt_fits(img, norm = LogNorm()):
    plt.imshow(img, norm=LogNorm(),origin='low')   
    plt.colorbar()
    plt.show()
    
