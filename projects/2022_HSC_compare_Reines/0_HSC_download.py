#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 15:27:48 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import glob

#%%Cutout using online cutout tool:
from galight.hsc_utils import hsc_image, hsc_psf
import os

#Read files:
f = open("2021_previous/Reines_ID_RaDec.txt","r")
string = f.read()
# object_id,ra,dec='150200.03+444541.9', 225.5694427, 2.955502
lines = string.split('\n')   # Split in to \n

dr='dr3'
# rerun='s21a_dud'
rerun='s19a_wide'

bands = 'GRIZY'  #Band that will be download

for i in range(len(lines)):
    object_id,ra,dec=lines[i].split(' ')
    ra,dec = float(ra), float(dec) 
    print('Downloading data with PSF... ... ...', i)
    out_dir='./s19a/' + object_id
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    try:
        hsc_image.get_cutouts(object_id,ra,dec,out_dir,dr=dr,rerun=rerun,filters=bands,fov_arcsec=120)
        hsc_psf.get_psfs(object_id,ra,dec,out_dir,dr=dr,rerun=rerun,filters=bands)
    except:
        None
    files = glob.glob(out_dir+'/*fits')
    if files == []:
        os.rmdir(out_dir)