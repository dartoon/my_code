#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 17:54:51 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import numpy as np
line_means = ['id', 'z', 'mbh', 'mbh_err', 'mgal', 'ps_mag', 'spectra', 'bit']
infers  = np.loadtxt('HSC_fitting/sdss_quasar_mbh.txt', dtype=str)
IDs_ = infers[:, 0]
HSC_z_ = infers[:,1].astype(np.float)
HSC_Mstar_ = infers[:,4].astype(np.float)
HSC_MBHs_ = infers[:,2].astype(np.float)
HSC_MBHs_err_ = infers[:,3].astype(np.float)
HSC_label_ = infers[:,-1]
ps_mag_ = infers[:,5].astype(np.float)

HSC_z, HSC_Mstar, HSC_MBHs, ps_mag = [], [], [], []
for i in range(len(IDs_)):
    if HSC_label_[i] in ['eboss_core', 'boss_core', 'ugri']:
        HSC_z.append(HSC_z_[i])
        HSC_Mstar.append(HSC_Mstar_[i])
        HSC_MBHs.append(HSC_MBHs_[i])
        ps_mag.append(ps_mag_[i])

# HSC_RA_ = infers[:,2].astype(np.float)
# HSC_DEC_ = infers[:,3].astype(np.float)