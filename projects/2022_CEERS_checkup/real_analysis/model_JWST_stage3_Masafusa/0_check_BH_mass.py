#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 21:30:00 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

f = open("target_idx_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n

idx_ = 0
line = lines[idx_]
idx, target_id, RA, Dec = line.split(' ')
RA, Dec = float(RA), float(Dec)

cata_file = '../../catalog_regions/SDSS_DR16Q_v4.fits'

# cata_file = '../../catalog_regions/spiders_quasar_bhmass-DR16-v1.fits'
hdul = pyfits.open(cata_file)
table = hdul[1].data
name_ = hdul[1].columns


# #%%
# for j in range(len(table)):
#     name_, RA_, Dec_ = table[j][1], table[j][2], table[j][3]
#     if np.sqrt((RA_- RA)**2 + (Dec_ - Dec)**2)*3600 < 100:
#         print(name_)