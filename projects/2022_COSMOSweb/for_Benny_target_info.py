#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:43:10 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import warnings, pickle
warnings.filterwarnings("ignore")
filt_i = 1
filt = ['F115W', 'F150W','F277W', 'F444W'][filt_i]
cata_list = pickle.load(open('material/cata_list.pkl','rb'))

report_ids = ['cid_1188', 'cid_1203', 'cid_1240',  'cid_1522', 'cid_1551', 'cid_2470']

for tid in report_ids:
    for clist in cata_list:
        if tid == clist[-1]:
            print(tid, clist[1], clist[2])