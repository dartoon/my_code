#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 10:13:22 2020

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

#%%This is the result of the SED fitting parameters
filename_p  = 'summary_20210719.fits'
filename_p = glob.glob(filename_p)[0]
hdul = pyfits.open(filename_p)
table = hdul[1].data
name = hdul[1].columns
age_idx = [i for i in range(len(name)) if 'AGE0' in str(name[i])][0]
age = str(round(10**table[1][age_idx],2)) # 'AGE:' Gyr
# z_idx = [i for i in range(len(name)) if 'zmc' in str(name[i])][0]
# z = str(round(table[1][z_idx],2))
mel_idx = [i for i in range(len(name)) if 'Z0' in str(name[i])][0]
mel = str(round(10**table[1][mel_idx],2)) # 'Mel:' Z*/Z_sun

filename_p  = 'SFH_20210719.fits'
filename_p = glob.glob(filename_p)[0]
hdul = pyfits.open(filename_p)
table1 = hdul[1].data
name1 = hdul[1].columns
age_idx = [i for i in range(len(name1)) if 'time' in str(name1[i])][0]
age = str(round(table1[0][age_idx],2)) # 'AGE:' Gyr

Mstel_id = [i for i in range(len(name1)) if 'Mstel50' in str(name1[i])][0]
Mstel = str(round(table1[0][Mstel_id],2)) # 'Mel:' Z*/Z_sun

# file = 'SFH_20210719.fits'
# fd_sfh = pyfits.open(file)[1].data

# file_sum = 'summary_20210719.fits'
# hd_sum = pyfits.open(file_sum)[1].header
# fd_sum = pyfits.open(file_sum)[1].data
# #print(hd_sum)

# time = fd_sfh['time']
# sfr16 = fd_sfh['SFR16']
# sfr50 = fd_sfh['SFR50']
# sfr84 = fd_sfh['SFR84']
# #sfr16 = fd_sfh['SFR16']
# #sfr50 = fd_sfh['SFR50']
# #sfr84 = fd_sfh['SFR84']
# # Plot SFH;
# plt.errorbar(time, sfr50, yerr=[sfr50-sfr16,sfr84-sfr50])

# plt.xscale('log')
# plt.xlabel('Lookback time / Gyr')
# plt.ylabel('SFR / Msun yr$^{-1}$')



#%%This is the spectra template
spec_file = "./gsf_spec_*.fits"
spec_file = glob.glob(spec_file)[0]
spec1d = pyfits.open(spec_file)  
table_spec = spec1d[1].data
plt.figure(figsize=(10, 6))
f = table_spec['f_model_50']            #This is the best-fit spectra (i.e. 50% by MCMC).
wave = table_spec['wave_model']/10000.
plt.plot(wave, f)

f = table_spec['f_model_16']
wave = table_spec['wave_model']/10000.
plt.plot(wave, f)

f = table_spec['f_model_84']
wave = table_spec['wave_model']/10000.
plt.plot(wave, f)

#Load the filter information:
f200w_fil = np.loadtxt('../354.fil')
f200w_fil[:,2] = f200w_fil[:,2]/f200w_fil[:,2].max()* f.max()/6
plt.plot(f200w_fil[:,1]/10000., f200w_fil[:,2], label='JWST F200W filter response', c='blue')

f356w_fil = np.loadtxt('../356.fil')
f356w_fil[:,2] = f356w_fil[:,2]/f356w_fil[:,2].max()* f.max()/6
plt.plot(f356w_fil[:,1]/10000., f356w_fil[:,2], label='JWST F356W filter response', c='orange')

#The host magnitude infered from 2D decomposition
mag = [24.4027, 23.0435]  #The input AB magnitude in F200W and F356W

fnu = [10 ** ((mag[i]-25)/(-2.5)) for i in range(len(mag))]
lam = np.array([np.mean(f200w_fil[:,1]), np.mean(f356w_fil[:,1])])

mag = -2.5*np.log10(fnu) + 25
flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**19
lam = lam/10000.
plt.scatter(lam, flambda, c='r', s=100, zorder =100)

plt.xlim(0.25, 8)
xmin, xmax, ymin, ymax = plt.axis()
plt.title('SED inference for J0859+0022', fontsize=25)
plt.text( (xmax-xmin)*0.5, (ymax-ymin)*0.6, 'inferred stellar population:', fontsize=17)
plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.53, 'age={0} Gyr'.format(age), fontsize=17)
plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.46, 'Z$_*$/Z$_\odot$={0}'.format(mel), fontsize=17)
plt.text( (xmax-xmin)*0.5, (ymax-ymin)*0.39, 'inferred stellar mass:', fontsize=17)
plt.text( (xmax-xmin)*0.77, (ymax-ymin)*0.32, 'M$_*$/M$_\odot$=10^'+Mstel, fontsize=17)
plt.legend(prop={'size':15}, loc = 1)
plt.tick_params(labelsize=15)
plt.xlabel("um",fontsize=27)
plt.ylabel(r"f$_\lambda$ 10$^{\rm -19}$ erg/s/cm$^2$/$\AA$",fontsize=27)
#plt.yticks([])
# plt.savefig('host_stellar_sed.pdf')
plt.show()

#%%
#Calculate the AB mag based on the SED and the filter:
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

#Calculate AB mag based on this spec and filter:
array_spec = np.array([list(wave*10000), list(f*10**(-19))]).T
f200w_filt_flam = cal_filt_flam(array_spec, f200w_fil[:,1:])
f200w_lam = np.mean(f200w_fil[:,1])
f200w_mag = -2.5 * np.log10(f200w_filt_flam) - 2.402 - 5.0 * np.log10(f200w_lam)
print("F200w filter AG magnitude:", f200w_mag)

array_spec = np.array([list(wave*10000), list(f*10**(-19))]).T
f356w_filt_flam = cal_filt_flam(array_spec, f356w_fil[:,1:])
f356w_lam = np.mean(f356w_fil[:,1])
f356w_mag = -2.5 * np.log10(f356w_filt_flam) - 2.402 - 5.0 * np.log10(f356w_lam)
print("F356w filter AG magnitude:", f356w_mag)

