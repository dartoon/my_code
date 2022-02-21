#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 10:42:19 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle
import glob

target_info = pickle.load(open('target_info.pkl','rb'))
res_files = glob.glob('sim_results/*filt_f35*pkl')
res_files.sort()
for file in res_files:
    res = pickle.load(open(file,'rb'))
    # print(round(res['true_host_flux'],1), round(res['inferred_host_flux'],1))
    print(res['PSF_id_true'], res['PSF_id_model'])
    # print(res['host_Reff_kpc'], res['host_flux_ratio'])

    