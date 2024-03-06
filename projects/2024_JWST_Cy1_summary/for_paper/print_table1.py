#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:07:58 2024

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info

for idx in range(10):
    info = target_info[str(idx)]
    target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
    import glob
    folder = '../../2022_JWST_QSOz6/NIRCam_data/*/bkg_removed/'   #!!!
    jwst_all_filenames = glob.glob(folder+'*{0}*{1}*{2}*_rmbkg.fits'.format(target_id[:5],target_id[-4:], 'F356W'))  #For NIRCam
    jwst_all_filenames.sort()
    file = jwst_all_filenames[0]
    fitsfile = pyfits.open(file)
    #Target Name, RA, Dec, redshift
    # print(target_id, RA, Dec, z, fitsfile[0].header['DATE-OBS'], fitsfile[0].header['EFFEXPTM'])
    target_name = fitsfile[0].header['TARGNAME']
    target_name = target_name.replace('-', '$-$')
    
    TARG_RA = '{0:.5f}'.format(fitsfile[0].header['TARG_RA']) #round(fitsfile[0].header['TARG_RA'],5)
    
    TARG_DEC = '{0:.5f}'.format(fitsfile[0].header['TARG_DEC'])  #str(round(fitsfile[0].header['TARG_DEC'],5))
    TARG_DEC = TARG_DEC.replace('-', '$-$')
    print(target_name, '&', TARG_RA, '&', 
          TARG_DEC, '&', z, '&', fitsfile[0].header['DATE-OBS'], '&', 
          fitsfile[0].header['EFFEXPTM'], '\\\\')
