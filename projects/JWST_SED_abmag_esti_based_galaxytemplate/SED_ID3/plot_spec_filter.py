#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 10:13:22 2020

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
# from specutils import Spectrum1D
import glob

filename_p  = 'SFH_0859_PA00_param.fits'
hdul = pyfits.open(filename_p)
table = hdul[1].data
name = hdul[1].columns
age_idx = [i for i in range(len(name)) if 'T_MW' in str(name[i])][0]
age = str(round(10**table[1][age_idx],3)) # 'AGE:' Gyr
z_idx = [i for i in range(len(name)) if 'zmc' in str(name[i])][0]
z = str(round(table[1][z_idx],2))
mel_idx = [i for i in range(len(name)) if 'Z_MW' in str(name[i])][0]
mel = str(round(10**table[1][mel_idx],3)) # 'Mel:' Z*/Z_sun

Mstel_id = [i for i in range(len(name)) if 'Mstel' in str(name[i])][0]
Mstel = table[1][Mstel_id] # The scaled stellar mass in the template

SFR_id = [i for i in range(len(name)) if 'SFR' in str(name[i])][0]
SFR = str(round(table[1][SFR_id],2)) # unit = 'Msun/yr'

M_dyn = 56. #For ID 3
s_mass = np.log10(M_dyn * 10 **9) 
ratio = 10**s_mass/ 10**Mstel

fnu_ini = 0.5 #initial fnu set as 0.5
fnu = fnu_ini * ratio
zp_fnu = 25  # zp_fnu == 25
ab_mag = -2.5*np.log10(fnu) + zp_fnu 
print("HSC-Y filter AG magnitude:", ab_mag)

#%%
spec1d = pyfits.open("./gsf_spec_0859.fits")  
table_spec = spec1d[1].data
plt.figure(figsize=(10, 6))
f = table_spec['f_model_50']
f = f * ratio # re-scale the f
wave = table_spec['wave_model']/10000.

#%%
plt.figure(figsize=(8, 6))
plt.plot(wave, f)

hsc_y_fil = np.loadtxt('hsc_y.fil')
hsc_y_fil[:,2] = hsc_y_fil[:,2]/hsc_y_fil[:,2].max()* f.max()/6
plt.plot(hsc_y_fil[:,1]/10000., hsc_y_fil[:,2], label='HSC-Y response')

plt.xlabel("um",fontsize=27)
plt.ylabel(r"f$_\lambda$ 10$^{\rm -19}$ erg/s/cm$^2$/$\AA$",fontsize=27)
plt.tick_params(labelsize=12)
plt.legend(prop={'size':15})
plt.xscale("log")
plt.xlim(0.5, 20 )
xmin, xmax, ymin, ymax = plt.axis()

plt.text( (xmax-xmin)*0.1, (ymax-ymin)*0.75, 'stellar population with:', fontsize=17)
plt.text( (xmax-xmin)*0.2, (ymax-ymin)*0.68, 'z={0}'.format(z), fontsize=17)
plt.text( (xmax-xmin)*0.2, (ymax-ymin)*0.6, 'age={0} Gyr'.format(age), fontsize=17)
plt.text( (xmax-xmin)*0.2, (ymax-ymin)*0.52, 'Z*/Z_sun={0} '.format(mel), fontsize=17)
plt.show()

#%%
#Calculate the Y mag based on the SED:
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

y_fil = np.loadtxt('hsc_y.fil')
y_fil = y_fil[:,1:]
array_spec = np.array([list(wave*10000), list(f*10**(-19))]).T
y_filt_flam = cal_filt_flam(array_spec, y_fil)
HSC_y_lam = 9775.1
y_mag = -2.5 * np.log10(y_filt_flam ) - 2.402 - 5.0 * np.log10(HSC_y_lam)
print("HSC-Y filter AG magnitude:", y_mag)
