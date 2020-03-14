#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 10:13:22 2020

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from specutils import Spectrum1D
import glob

#%%To estimate the AB mag for a given stellar template
filename_p  = 'Spectra_at_z1/SFH_*_PA00_param.fits'
filename_p = glob.glob(filename_p)[0]
hdul = pyfits.open(filename_p)
table = hdul[1].data
name = hdul[1].columns
age_idx = [i for i in range(len(name)) if 'T_MW' in str(name[i])][0]
age = str(round(10**table[1][age_idx],3)) # 'AGE:' Gyr
z_idx = [i for i in range(len(name)) if 'zmc' in str(name[i])][0]
z = str(round(table[1][z_idx],2))
mel_idx = [i for i in range(len(name)) if 'Z_MW' in str(name[i])][0]
mel = str(round(10**table[1][mel_idx],3)) # 'Mel:' Z*/Z_sun

spec1d  = pyfits.open("Spectra_at_z1/gsf_spec_101.fits")  
name_spec = spec1d[1].columns
table_spec = spec1d[1].data
array_spec = np.zeros((len(table_spec), 7))
for i in range(len(table_spec)):
    array_spec[i, :] = table_spec[i]

plt.figure(figsize=(10, 6))
array_spec[:,2] =  array_spec[:,2]/ array_spec[:,2].max() * 3.0
plt.plot(array_spec[:,0]/10000., array_spec[:,2])

#%%To estimate the AB mag for a given stellar template
filename_p  = 'Spectra_at_z2/SFH_*_PA00_param.fits'
filename_p = glob.glob(filename_p)[0]
hdul = pyfits.open(filename_p)
table = hdul[1].data
name = hdul[1].columns
age_idx = [i for i in range(len(name)) if 'T_MW' in str(name[i])][0]
age = str(round(10**table[1][age_idx],3)) # 'AGE:' Gyr
z_idx = [i for i in range(len(name)) if 'zmc' in str(name[i])][0]
z = str(round(table[1][z_idx],2))
mel_idx = [i for i in range(len(name)) if 'Z_MW' in str(name[i])][0]
mel = str(round(10**table[1][mel_idx],3)) # 'Mel:' Z*/Z_sun

spec1d  = pyfits.open("Spectra_at_z2/gsf_spec_102.fits")  
name_spec = spec1d[1].columns
table_spec = spec1d[1].data
array_spec = np.zeros((len(table_spec), 7))
for i in range(len(table_spec)):
    array_spec[i, :] = table_spec[i]

array_spec[:,2] =  array_spec[:,2]/ array_spec[:,2].max() * 3.0
plt.plot(array_spec[:,0]/10000., array_spec[:,2])
