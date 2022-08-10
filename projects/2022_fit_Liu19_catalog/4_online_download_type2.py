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
f = open("catalog/SDSS_DR17_AGN2_SNR2_z001035_objid_r90_catalog.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
import glob

no_per_core = int(len(lines)/4)

core = 3
for i, line in enumerate(lines[no_per_core*core+1:no_per_core*(core+1)+1]): 
    ID, Ra, Dec, z, petro = line.split(' ')
    out_dir = 'type2_online_data_download/' + ID
    if glob.glob('type2_online_data_download/'+ID) == []:
        os.makedirs(out_dir)
    if glob.glob('type2_online_data_download/'+ID+'/*fits') == []:
        hsc_image.get_cutouts(ID,Ra,Dec,out_dir,dr=dr,rerun=rerun,filters=bands,fov_arcsec=40)
        hsc_psf.get_psfs(ID,Ra,Dec,out_dir,dr=dr,rerun=rerun,filters=bands)
        print("core", core, i,ID, 'finished')
    else:
        print(i,ID, 'skip')
        