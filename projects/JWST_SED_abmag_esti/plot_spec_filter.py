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
#filename  = 'SED_files/summary_*_PA00.fits'
#filename = glob.glob(filename)[0]
#hdul_sum = pyfits.open(filename)
#name_sum = hdul_sum[1].columns
#table_sum = hdul_sum[1].data
#z_idx = [i for i in range(len(name_sum)) if 'zmc' in str(name_sum[i])][0]
#mel_idx = [i for i in range(len(name_sum)) if 'logZsun' in str(name_sum[i])][0]
#z = table_sum[1][z_idx] #Redshift
#mel = table_sum[1][mel_idx] #Metallicity 
filename_p  = 'SED_files/SFH_*_PA00_param.fits'
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

#%%
spec1d = hdul_sum = pyfits.open("gsf_spec_1.fits")  
name_spec = spec1d[1].columns
table_spec = spec1d[1].data
array_spec = np.zeros((len(table_spec), 7))
for i in range(len(table_spec)):
    array_spec[i, :] = table_spec[i]

#%%

sov_jwst_f144w_fil = np.loadtxt('JWST_NIRCam.F444W.dat')
sov_jwst_f144w_fil[:,1] = sov_jwst_f144w_fil[:,1]/sov_jwst_f144w_fil[:,1].max()* array_spec[:,2].max()/6

file_text = []
for i in range(len(sov_jwst_f144w_fil)):
	file_text.append([i+1, sov_jwst_f144w_fil[i,0], sov_jwst_f144w_fil[i,1]])
file_text = np.asarray(file_text)

plt.figure(figsize=(8, 6))
plt.plot(array_spec[:,0]/10000., array_spec[:,2])
plt.plot(sov_jwst_f144w_fil[:,0]/10000., sov_jwst_f144w_fil[:,1], label='NIRCam F444W response curve')
plt.xlim(0.25, 8)
xmin, xmax, ymin, ymax = plt.axis()
plt.text( (xmax-xmin)*0.45, (ymax-ymin)*0.75, 'stellar population with:', fontsize=17)
plt.text( (xmax-xmin)*0.7, (ymax-ymin)*0.68, 'z={0}'.format(z), fontsize=17)
plt.text( (xmax-xmin)*0.7, (ymax-ymin)*0.6, 'age={0} Gyr'.format(age), fontsize=17)
plt.text( (xmax-xmin)*0.7, (ymax-ymin)*0.52, 'Z*/Z_sun={0} '.format(mel), fontsize=17)
plt.legend(prop={'size':15})
plt.tick_params(labelsize=12)
plt.xlabel("um",fontsize=27)
#plt.yticks([])