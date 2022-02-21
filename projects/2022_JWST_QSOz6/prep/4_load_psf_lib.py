#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 17:22:06 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle
psfs, FWHMs = pickle.load(open('f150w_psfs.pkl','rb'))
# psfs, FWHMs = pickle.load(open('f356w_psfs.pkl','rb'))

from galight.tools.astro_tools import plt_many_fits
for i in range(len(psfs)):
    # plt_many_fits(psfs[i], FWHMs[i], 'FWHM')
    plt_many_fits(psfs[i], [np.sum(psfs[i][j]) for j in range(len(psfs[i]))], 'Flux')
#     if len(psfs[i]) != len(FWHMs[i]):   
#           print(i, "Shape not right for!")