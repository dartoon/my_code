#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 16:50:05 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle


filt = 'f356w'
_psfs, _FWHMs = pickle.load(open(filt+'_psfs.pkl','rb'))
psfs, FWHMs, fluxs = [], [], []
for i in range(len(_psfs)):
    psfs = psfs + _psfs[i]
    FWHMs = FWHMs + _FWHMs[i]
    fluxs = fluxs + [np.sum(_psfs[i][j]) for j in range(len(_psfs[i]))]
    
FWHMs, fluxs = np.array(FWHMs), np.array(fluxs)

def find_close_PSF_idx(psf_id):
    sort = np.argsort(abs(FWHMs[psf_id] - FWHMs))[1:]
    idx = sort[0] #Setting a initial value
    for i in sort:
        # print( abs((fluxs[i]- fluxs[psf_id])/fluxs[psf_id]  ) )
        if abs( (fluxs[i]- fluxs[psf_id])/fluxs[psf_id])<0.5 and fluxs[i]>500:
            idx = i
            break
    # print(sort)
    return idx
# print(find_close_PSF_idx(0))