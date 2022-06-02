#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 10:06:11 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

#Select which data you plan to download
dr='dr4'
# rerun='s21a_dud'  #Deep or UltraDeep
rerun='s21a_wide'  #Wide 
#Select how many bands you plan to Download (G R I Z Y)
bands = 'GRIZY'  #Band that will be download
import os
f = open("catalog.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
import glob
from galight.hsc_utils import hsc_image, hsc_psf
dr='dr4'
rerun='s21a_wide'  #Wide 
ct = 0
from galight.data_process import DataProcess

for _, line in enumerate(lines):
    ID, Ra, Dec = line.split(' ')
    for b in bands:
        img_globname = glob.glob('/Volumes/Seagate_Expansion_Drive/data_backup/Liu19_catalog/gfarm_data_download/{0}_HSC-{1}.fits'.format(ID,b)) + glob.glob(
            '/Volumes/Seagate_Expansion_Drive/data_backup/Liu19_catalog/online_data_download/{0}/*cutout*-{1}-*.fits'.format(ID,b) )
        psf_globname = glob.glob('/Volumes/Seagate_Expansion_Drive/data_backup/Liu19_catalog/gfarm_data_download/{0}_HSC-{1}_psf.fits'.format(ID,b)) + glob.glob(
            '/Volumes/Seagate_Expansion_Drive/data_backup/Liu19_catalog/online_data_download/{0}/*psf*-{1}-*.fits'.format(ID,b) )
        glob_files = glob.glob('fit_result/{0}*-{1}.pkl'.format(ID,b))
        if len(img_globname) >= 1 and len(psf_globname) >= 1:
            if glob_files == []:
                print(ID, b)
        # ['10547', '13120', '13338', '14503', '1956', '2831', '5328', '7007', '8457', '8917', '9009', '9209']