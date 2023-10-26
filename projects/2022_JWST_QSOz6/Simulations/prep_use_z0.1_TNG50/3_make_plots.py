#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 08:12:06 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle
import glob

res_files = glob.glob('results_extensive_sim/qsoID0_filt_f*pkl') + glob.glob('results_extensive_sim/qsoID1_filt_f*pkl')
# res_files = glob.glob('sim_results/*seed2_small*_result.pkl')
res_files.sort()

imgs = []
for file in res_files:
    res = pickle.load(open(file,'rb'))
    imgs.append(res.target_stamp)

