#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 17:21:06 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

cata_file = './TNG50_catalog/TNG50-1_091_Catalogue.csv'

import pandas as pd
cata = pd.read_csv(cata_file)

smass = cata['SubhaloMassType_stars']
flag = cata['Unnamed: 0']

reff_star = cata['SubhaloHalfmassRadType_stars']
reff_dm = cata['SubhaloHalfmassRadType_dm']

# print(flag[(smass>10.4)*(smass<10.5)*(reff_star<3)*(reff_star>1)])

bool_ = (smass>10.3)*(smass<10.8)*(reff_star<2)*(reff_star>0.5)
# print(flag[bool_])

for k in flag[bool_]:
    print('http://idark.ipmu.jp/hsc405/SKIRT9/Photometry/091/shalo_091-{0}_v0_photo.fits'.format(k))
