#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 18:08:19 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob

all_IDs  = glob.glob('../_John_fitted/*/proof-BHBH.pdf')
all_IDs = [all_IDs[i].split('/')[2].split('_HSC')[0] for i in range(len(all_IDs))]

fitted_IDs = glob.glob('../proofBHBH/model_Iband_zover_1/*')
fitted_IDs = [fitted_IDs[i].split('/')[-1] for i in range(len(fitted_IDs))]

from tools import read_info
from astropy.coordinates import SkyCoord
from astropy import units as u
write_file = open('proof-BHBH.txt','w') 
for ID in all_IDs:
    # if ID not in fitted_IDs:
    RA, Dec, z = read_info(ID)
    if RA == -99:
        pos = SkyCoord('{0} {1}'.format(ID[:2]+':'+ID[2:4]+':'+ID[4:9], ID[9:12]+':'+ID[12:14]+':'+ID[14:]), unit=(u.hourangle, u.deg))
        RA, Dec = pos.ra.degree, pos.dec.degree
    write_file.write(ID+" {0} {1} {2}\n".format(RA, Dec, z))
write_file.close()