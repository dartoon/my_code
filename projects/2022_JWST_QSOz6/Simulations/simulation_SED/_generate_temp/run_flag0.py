#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 09:39:13 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#%%
from gsf import gsf
gsf.run_gsf_all('sample.input', 1, idman=None)
plt.show()

#%%
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
import glob
spec1d_file = glob.glob('./gsf_spec_*.fits')[0]
spec1d = pyfits.open(spec1d_file)  
name_spec = spec1d[1].columns
unit =  spec1d[1].header['TUNIT3']
table_spec = spec1d[1].data

steller_file = glob.glob('./SFH_*.fits')[0]
hdul = pyfits.open(steller_file)
info1 = hdul[0].header 
print('redshift', float(info1['ZMC_50']))
print('smass', info1['Mstel_50'], info1['Mstel_16'], info1['Mstel_84']) 

wave = table_spec['wave_model']
f_50 = table_spec['f_model_50']
plt.plot(wave/10000., f_50, 'black', alpha=0.7)
plt.xlim(0.25, 6)
plt.ylim(0., 0.88)

#%%
jwst_filt_id = {'F115W': '352', 'F150W': '353', 'F200W': '354', 
            'F277W': '355', 'F356W': '356', 'F444W': '357', 'F410M': '362'}
filt = 'F356W'
fid = jwst_filt_id[filt]
f_fil = np.loadtxt('../../../../template/gsf_temp/filter/{0}.fil'.format(fid))
lam = np.median(f_fil[1:,1])
f_array = np.vstack((wave, f_50)).T
filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
flambda = filt_flam 
# flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**float(unit.split('e-')[1][:2])
mag= -2.5*np.log10(flambda / (10**float(unit.split('e-')[1][:2]))) - (2.402+5.0*np.log10(lam)) 

