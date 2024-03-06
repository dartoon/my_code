#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 23:46:56 2024

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info
from galight.tools.astro_tools import read_pixel_scale

filt = 'F150W' #!!!

for i in range(len(target_info)):
    idx = i #!!! 
    info = target_info[str(idx)]
    target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
    
    import glob
    folder = '../../2022_JWST_QSOz6/NIRCam_data/*/bkg_removed/'   #!!!
    jwst_all_filenames = glob.glob(folder+'*{0}*{1}*{2}*_rmbkg.fits'.format(target_id[:5],target_id[-4:], filt))  #For NIRCam
    jwst_all_filenames.sort()
    print(i)
    print(jwst_all_filenames)
    
    file = jwst_all_filenames[0]
    _fitsFile = pyfits.open(file)
    # fov_image = _fitsFile[1].data # check the back grounp
    header = _fitsFile[1].header
    print(read_pixel_scale(header))