#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 22:59:29 2022

@author: Dartoon

For Benny Trakhtenbrot
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'..')
from def_functions import RA_Dec_in_fit
import warnings
warnings.filterwarnings("ignore")

RA_DEC_list = [[150.004380, 2.038898],
[149.869660, 2.294046],
[150.208850, 2.481935],
[150.297250, 2.148846],
[149.472888, 2.793379],
[150.735585, 2.199557],
[150.240801, 2.659037],
[150.737172, 2.722557],
[150.620069, 2.671382],
[149.529103, 2.380143],
[150.720703, 2.693635],
[150.767390, 2.739021],
[150.705563, 2.629612],
[150.058920, 2.015179]]

import glob

folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_JWST_Masafusa'
filenames = glob.glob(folder+'/bkg_removed/'+'*.fits')

for RA, Dec in RA_DEC_list:
    print(RA, Dec)
    files = RA_Dec_in_fit(all_files=filenames, RA=float(RA), Dec=float(Dec))
    print(files)
    