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

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

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
filename_p  = './SFH_*_PA00_param.fits'
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

Mstel_id = [i for i in range(len(name)) if 'Mstel' in str(name[i])][0]
Mstel = str(round(table[1][Mstel_id],2)) # 'Mel:' Z*/Z_sun

#%%
spec1d = pyfits.open("./gsf_spec_1234.fits")  
table_spec = spec1d[1].data
plt.figure(figsize=(10, 6))
f = table_spec['f_model_50']
wave = table_spec['wave_model']/10000.
plt.plot(wave, f)

f = table_spec['f_model_16']
wave = table_spec['wave_model']/10000.
plt.plot(wave, f)

f = table_spec['f_model_84']
wave = table_spec['wave_model']/10000.
plt.plot(wave, f)

# plt.plot(array_spec[:,0]/10000., array_spec[:,2])

f200w_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf-version1.4/example/filter/354.fil')
f200w_fil[:,2] = f200w_fil[:,2]/f200w_fil[:,2].max()* f.max()/6
plt.plot(f200w_fil[:,1]/10000., f200w_fil[:,2], label='JWST F200W filter response', c='blue')

f356w_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf-version1.4/example/filter/356.fil')
f356w_fil[:,2] = f356w_fil[:,2]/f356w_fil[:,2].max()* f.max()/6
plt.plot(f356w_fil[:,1]/10000., f356w_fil[:,2], label='JWST F356W filter response', c='orange')


fnu = np.array([1.733, 6.062])
lam = np.array([np.mean(f200w_fil[60:,1]), np.mean(f356w_fil[60:,1])])
yerr_l = np.array([0.03, 0.03])

mag = -2.5*np.log10(fnu) + 25
flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**19
lam = lam/10000.
plt.scatter(lam, flambda, c='r', zorder =100)
plt.errorbar(lam, flambda, yerr=yerr_l*flambda,fmt='.',color='gray',markersize=1, zorder =90)

#plt.plot(sov_jwst_f144w_fil[:,0]/10000., sov_jwst_f144w_fil[:,1], label='NIRCam F444W response curve')
plt.xlim(0.25, 8)
xmin, xmax, ymin, ymax = plt.axis()
plt.title('SED inference for J0859+0022', fontsize=25)
# plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.6, 'z={0}'.format(z), fontsize=17)
plt.text( (xmax-xmin)*0.5, (ymax-ymin)*0.6, 'inferred stellar population:', fontsize=17)
plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.53, 'age={0:.3f} Gyr'.format(float(age)), fontsize=17)
plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.46, 'Z$_*$/Z$_\odot$={0}'.format(mel), fontsize=17)
plt.text( (xmax-xmin)*0.5, (ymax-ymin)*0.39, 'inferred stellar mass:', fontsize=17)
plt.text( (xmax-xmin)*0.77, (ymax-ymin)*0.32, 'M$_*$/M$_\odot$=10^'+Mstel, fontsize=17)
plt.legend(prop={'size':15}, loc = 1)
plt.tick_params(labelsize=15)
plt.xlabel("um",fontsize=27)
plt.ylabel(r"f$_\lambda$ 10$^{\rm -19}$ erg/s/cm$^2$/$\AA$",fontsize=27)
#plt.yticks([])
plt.savefig('host_stellar_sed.pdf')