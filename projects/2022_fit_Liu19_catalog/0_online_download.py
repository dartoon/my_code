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
f = open("Liu19_HSC_DR4_objid_catalog.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
import glob
# for line in lines[:500]:
# for line in lines[500:1000]:
# for line in lines[1000:1500]:
for line in lines[1500:]:
    ID, Ra, Dec = line.split(' ')
    if glob.glob('gfarm_data_download/'+ID+'_HSC*fits') == []:
        print("no exit",ID)
        out_dir = './' + ID
        os.makedirs(ID)
        hsc_image.get_cutouts(ID,Ra,Dec,out_dir,dr=dr,rerun=rerun,filters=bands,fov_arcsec=120)
        hsc_psf.get_psfs(ID,Ra,Dec,out_dir,dr=dr,rerun=rerun,filters=bands)
    else:
        print('exit',ID)
        