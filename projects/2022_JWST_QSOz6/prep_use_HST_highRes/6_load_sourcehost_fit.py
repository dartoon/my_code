#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 21:43:01 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob,pickle

files = glob.glob('test_int/*seed1*pkl')
files.sort()
for file in files:
    res = pickle.load(open(file,'rb'))
    print(res['source_id'], res['inferred_n_sersic_same_psf'])