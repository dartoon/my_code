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
res_files = glob.glob('/Users/Dartoon/Downloads/sim_results/qsoID*filt_f35*pkl')
res_files.sort()

res_scatter_diffPSF, res_scatter_samePSF = [], []

ratio = []
for file in res_files:
    res = pickle.load(open(file,'rb'))
    # print(round(res['true_host_flux'],1), round(res['inferred_host_flux'],1))
    # print(res['PSF_id_true'], res['PSF_id_model'])
    # print(res['host_Reff_kpc'], res['host_flux_ratio'])
    # res_scatter_diffPSF.append( res['inferred_magnitude_diff_psf'] - res['true_host_mag'] )
    # res_scatter_samePSF.append( res['inferred_magnitude_same_psf'] - res['true_host_mag'] )
    
    res_scatter_diffPSF.append( res['inferred_R_sersic_diff_psf'] / res['galfit_Re'] )
    res_scatter_samePSF.append( res['inferred_R_sersic_same_psf'] / res['galfit_Re'] )
    ratio.append( res['true_host_flux_ratio'] )

#%%
plt.figure(figsize=(11,11))
plt.scatter(ratio, res_scatter_diffPSF)
plt.scatter(ratio, res_scatter_samePSF)
plt.ylim([0, 3])
# plt.ylim([-1.5, 1.5])
plt.show()