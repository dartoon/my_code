#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 15:09:22 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from galight.hsc_utils import hsc_image, hsc_psf

#Select which data you plan to download
dr='dr4'
# rerun='s21a_dud'  #Deep or UltraDeep
rerun='s21a_wide'  #Wide 
#Select how many bands you plan to Download (G R I Z Y)
bands = 'GRIZY'  #Band that will be download
import os
f = open("catalog/GAMA_BLAGN_sample_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
import glob

# for i, line in enumerate(lines[1:100]):
# for i, line in enumerate(lines[100:200]):
# for i, line in enumerate(lines[200:300]):
for i, line in enumerate(lines[300:]):    
    ID, Ra, Dec = line.split(',')
    out_dir = 'GAMA_online_data_download/' + ID
    if glob.glob('GAMA_online_data_download/'+ID) == []:
        os.makedirs(out_dir)
    if glob.glob('GAMA_online_data_download/'+ID+'/*fits') == []:
        hsc_image.get_cutouts(ID,Ra,Dec,out_dir,dr=dr,rerun=rerun,filters=bands,fov_arcsec=40)
        hsc_psf.get_psfs(ID,Ra,Dec,out_dir,dr=dr,rerun=rerun,filters=bands)
        print(i,ID, 'finished')
    else:
        print(i,ID, 'skip')
        