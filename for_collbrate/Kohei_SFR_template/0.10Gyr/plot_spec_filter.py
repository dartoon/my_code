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
filename_p  = './summary_20210802*.fits'
filename_p = glob.glob(filename_p)[0]
hdul = pyfits.open(filename_p)
table = hdul[1].data
name = hdul[1].columns
z_idx = [i for i in range(len(name)) if 'z_cz' in str(name[i])][0]
z = str(round(table[1][z_idx],2))
mel_idx = [i for i in range(len(name)) if 'Z0' in str(name[i])][0]
mel = str(round(10**table[1][mel_idx],2)) # 'Mel:' Z*/Z_sun

#%%
filename_p  = './SFH_20210802*.fits'
filename_p = glob.glob(filename_p)[0]
hdul = pyfits.open(filename_p)
table = hdul[1].data
name = hdul[1].columns
age_idx = [i for i in range(len(name)) if 'time' in str(name[i])][0]
age = table[0][age_idx] # 'AGE:' Gyr
Mstel_id = [i for i in range(len(name)) if 'Mstel50' in str(name[i])][0]
Mstel = table[0][Mstel_id] # 'Mstel


Mstel_nor = 7
nor = 10**Mstel_nor/10**Mstel

SFR_id = [i for i in range(len(name)) if 'SFR50' in str(name[i])][0]
SFR = table[0][SFR_id] # 'Mstel
print( np.log10(10**SFR*nor) )

#%%This is the spectra template
spec_file = "./gsf_spec_*.fits"
spec_file = glob.glob(spec_file)[0]
spec1d = pyfits.open(spec_file)  
table_spec = spec1d[1].data
plt.figure(figsize=(10, 6))
f = table_spec['f_model_50'] * nor            #This is the best-fit spectra (i.e. 50% by MCMC).
wave = table_spec['wave_model']/10000.
plt.plot(wave, f)

#Load the filter information:
f200w_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf/example/filter/354.fil')
f200w_fil[:,2] = f200w_fil[:,2]/f200w_fil[:,2].max()* f.max()/6
# plt.plot(f200w_fil[:,1]/10000., f200w_fil[:,2], label='JWST F200W filter response', c='blue')

f356w_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf/example/filter/356.fil')
f356w_fil[:,2] = f356w_fil[:,2]/f356w_fil[:,2].max()* f.max()/6
plt.plot(f356w_fil[:,1]/10000., f356w_fil[:,2], label='JWST F356W filter response', c='orange')

#The host magnitude infered from 2D decomposition
mag = [24.4027, 23.0435]  #The input AB magnitude in F200W and F356W

fnu = [10 ** ((mag[i]-25)/(-2.5)) for i in range(len(mag))]
lam = np.array([np.mean(f200w_fil[:,1]), np.mean(f356w_fil[:,1])])

mag = -2.5*np.log10(fnu) + 25
flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**19
lam = lam/10000.
plt.scatter(lam[1], flambda[1]* nor, c='r', s=100, zorder =100)

plt.xlim(0.25, 8)
plt.ylim(0,0.0005)
xmin, xmax, ymin, ymax = plt.axis()
plt.title('Spectrum', fontsize=25)
plt.text( (xmax-xmin)*0.5, (ymax-ymin)*0.6, 'inferred stellar population:', fontsize=17)
plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.53, 'age={0:.3f} Gyr'.format(age), fontsize=17)
plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.46, 'Z$_*$/Z$_\odot$={0}'.format(mel), fontsize=17)
plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.39, 'redshift {0}'.format(z), fontsize=17)
# plt.text( (xmax-xmin)*0.5, (ymax-ymin)*0.35, 'inferred stellar mass:', fontsize=17)

plt.text( (xmax-xmin)*0.77, (ymax-ymin)*0.32, 'M$_*$/M$_\odot$=10^{0:.1f}'.format(np.log10(10**Mstel*nor)), fontsize=17)
plt.text( (xmax-xmin)*0.67, (ymax-ymin)*0.25, 'SFR=10^{0:.3f} Msun/yr'.format(np.log10(10**SFR*nor) ), fontsize=17)


plt.legend(prop={'size':15}, loc = 1)
plt.tick_params(labelsize=15)
plt.xlabel("um",fontsize=27)
plt.ylabel(r"f$_\lambda$ 10$^{\rm -19}$ erg/s/cm$^2$/$\AA$",fontsize=27)
#plt.yticks([])
# plt.savefig('host_stellar_sed.pdf')

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

plt.text( (xmax-xmin)*0.1, (ymax-ymin)*0.85, 'F356W filter mag {0:.2f}'.format(f356w_mag), fontsize=17)
plt.show()