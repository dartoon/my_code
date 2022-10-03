#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 21:12:35 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob, pickle

filters =  ['F115','F150','F200','F277','F356', 'F410','F444']

for filt in filters:
    PSF_lib_files = glob.glob('material/*'+filt[:-1]+'*_PSF_Library.pkl')[0]
    PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
    if filt != 'F410':
        s = 'W'
    else:
        s = 'M'
        
    filt = filt+s
    
    print(filt, 'WCS Positions RA, Dec')
    for item in PSF_RA_DEC_list:
        print("{0:.4f} {1:.4f}".format(item[0], item[1]))
