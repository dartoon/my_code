#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 16:37:35 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pandas as pd
from scipy.interpolate import interp1d
def cal_filt_flam(array_spec, fil):
    '''
    Calculate a filter f_lambda
    Parameter
    --------
    Return
    --------
        Filter flux
    '''
    fil[:,1] /= fil[:,1].max()
    f_filt = interp1d(fil[:,0], fil[:,1], kind='cubic')
    int_spec = array_spec[(array_spec[:,0]>fil[0,0]) * (array_spec[:,0]<fil[-1,0])]
    int_flux = 0
    int_filt = 0
    for i in range(len(int_spec)-1):
        int_flux += int_spec[i,1] * f_filt(int_spec[i,0]) * (int_spec[i+1,0]-int_spec[i,0])
        int_filt += f_filt(int_spec[i,0]) * (int_spec[i+1,0]-int_spec[i,0])
    filt_flam = int_flux/ int_filt   
    return filt_flam

#%%
sample = pd.read_csv('HSCJ2255p0251_lines_edit.csv')
fitidx = 0
keys = [ 'line_Hb_na', 'line_O3_4959_c', 'line_O3_5008_c', 'line_O3_4959_w', 'line_O3_5008_w',]

sample = pd.read_csv('HSCJ2236p0032_lines_edit.csv')
fitidx = 1
keys = ['line_O3_4959_na', 'line_O3_5008_na', 'line_O3_4959_w',
        'line_O3_5008_w']


keys = ['line_Hb_br','line_O3_4959_na', 'line_O3_5008_na', 'line_O3_4959_w',
        'line_O3_5008_w']
keys = ['flux']
flux = sample[keys[0]]

for key in keys[1:]:
    flux += sample[key]
    
import sys
sys.path.insert(0,'../model_z6_data_id0')

from target_info import target_info
info = target_info[str(fitidx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

wave = sample['wavelength'] * (1+z)
flux = flux/ (1+z)
plt.figure(figsize=(10, 6))
plt.plot(wave, flux, label='Hbeta+OIII')
# plt.xlim([4600,5200])
hst_filt_id = {'F606W': '4', 'F814W':'6', 'F105W':'202', 'F125W':'203', 'F140W':'204', 'F160W':'205'}
jwst_filt_id = {'F115W': '352', 'F150W': '353', 'F200W': '354', 
            'F277W': '355', 'F356W': '356', 'F444W': '357', 'F410M': '362'}
filt_id = hst_filt_id | jwst_filt_id
f_fil = np.loadtxt('../../../template/gsf_temp/filter/{0}.fil'.format(filt_id['F356W']))        
f_array = np.vstack((wave, flux)).T
lam = np.median(f_fil[1:,1])
# lam = lam. 
filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
plt.scatter(lam, filt_flam, marker="d", zorder =90,  s=280, facecolors='none', edgecolors='c', linewidths=2, label='filter flux density')

xmin, xmax, ymin, ymax = plt.axis()
top = flux.max()
f_fil[:,2] = f_fil[:,2]/f_fil[:,2].max() * (ymax-ymin) * 0.1
plt.plot(f_fil[1:,1], f_fil[1:,2], label='F356W filter response')
plt.xlabel(r"$\lambda$ ($\AA$)",fontsize=25)
plt.ylabel(r"f$_\lambda$  (10$^{\rm" + " -{0}}}$".format(17)+" erg s$^{-1}$ cm$^{-2}$$\mathrm{\AA}^{-1}$)",fontsize=25)
plt.legend(prop={'size':18}, ncol=1, loc = 1)
plt.tick_params(labelsize=20)
plt.show()

mag =  -2.5 * np.log10( filt_flam * 10**(-17) ) - 2.402 - 5.0 * np.log10( lam )

import pickle
fit_file = '../model_z6_data_id0/stage3_all/fit_material/fit_run_idx0_F356W_CombPsfsNO_2_0.pkl'#+\
fit_run = pickle.load(open(fit_file,'rb'))
zp = fit_run.zp

NIRCam_flux = 10**(-0.4*(mag-zp))
print(NIRCam_flux) #This value is before correction.

#!!! For example, J2255 should be 6/248*155.  248 is the total in NIRSpec, 155 is total in NIRCam
#!!! For example, J2236 should be 1.97/309*234. 248 is the total in NIRSpec, 155 is total in NIRCam

