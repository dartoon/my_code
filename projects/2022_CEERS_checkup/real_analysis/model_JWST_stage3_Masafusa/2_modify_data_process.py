#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 15:10:43 2022

@author: Dartoon
"""


import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle
import copy

files = glob.glob('fit_material/data_process_idx*.pkl')
files.sort()
collect_info = []
i_list = []
for i in range(len(files)):
# for i in [2768]:
    _file = files[i]
    idx_info = _file.split('idx')[1].split('_')[0]
    filt_info = _file.split('_psf')[0].split('_')[-1]
    if int(idx_info) in [31, 32, 56]:
        print(i)
        i_list.append(i)
        # data_process = pickle.load(open(_file,'rb'))
        