#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 21:22:02 2020

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM

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

filename_p  = 'qso_template.fits' #Download from:http://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/the-agn-atlas
#filename_p  = 'seyfert1_template.fits'
#filename_p  = 'ngc1068_template.fits'
hdul = pyfits.open(filename_p)
table_spec = hdul[1].data
name = hdul[1].columns

array_spec = np.zeros((len(table_spec), len(table_spec[0])))
for i in range(len(table_spec)):
    array_spec[i, :] = table_spec[i]
    
array_spec_temp = copy.deepcopy(array_spec)
#%%Plot the Spec at the rest-frame
array_spec = copy.deepcopy(array_spec_temp)
plt.figure(figsize=(10, 6))    
plt.plot(array_spec[:,0], array_spec[:,1])
idx = np.where( abs(array_spec[:,0] - 1450.) == np.min(abs(array_spec[:,0] - 1450.)) )[0][0]
plt.scatter(1450, array_spec[idx,1], c='r')
plt.tick_params(labelsize=15)
plt.xlabel("A",fontsize=27)
plt.ylabel(r"f$_\lambda$ erg/s/cm$^2$/A",fontsize=27)
plt.close()

#%%Shift to redshift z = 6 and norm based on M1450 value
M1450, m_zAB, m_yAB, z = -23.82, 22.775, 22.942, 6.10
#M1450, m_zAB, m_yAB, z = -25.31, 21.827, 21.614, 6.37
#M1450, m_zAB, m_yAB, z = -24.09, 22.768, 23.649, 6.39
#M1450, m_zAB, m_yAB, z = -24.73, 22.121, 22.045, 6.2
#M1450, m_zAB, m_yAB, z = -24.69, 22.402, 22.326, 6.26

array_spec = copy.deepcopy(array_spec_temp)
array_spec[:,0] *= (1+z)
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dl = cosmo.luminosity_distance(z=z).value  
#m1450 = M1450 + 2.5*np.log10((dl*10**6/10)**2)   # Mv = m - 2.5 log[ (d/10)**2 ],  d in pc scale
m1450 = M1450 + 5*(np.log10(dl*10**6)-1)
lam = 1450 * (1+z)
f_lambda_1450 = 10**(-0.4*(m1450+2.402+5.0*np.log10(lam)))
idx = np.where( abs(array_spec[:,0] - lam) == np.min(abs(array_spec[:,0] - lam)) )[0][0]
norm = f_lambda_1450/array_spec[idx, 1]
array_spec[:,1] *= norm

z_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf/gsf/example/filter/317.fil')
z_fil = z_fil[:,1:]
z_fil[:,1] = z_fil[:,1]/z_fil[:,1].max()* array_spec[:,1].max()/6

y_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf/gsf/example/filter/318.fil')
y_fil = y_fil[:,1:]
y_fil[:,1] = y_fil[:,1]/y_fil[:,1].max()* array_spec[:,1].max()/6

plt.figure(figsize=(10, 6))    
plt.plot(array_spec[:,0], array_spec[:,1])
plt.plot(z_fil[:,0], z_fil[:,1], label='HSC Z filter response', c='green')
plt.plot(y_fil[:,0], y_fil[:,1], label='HSC Y filter response', c='orange')
plt.scatter(lam, f_lambda_1450, c='r', marker="*", s =200, zorder = 10, label='observed M1450')

HSC_z_lam = z_fil[np.where(z_fil[:,1] == z_fil[:,1].max())[0][0],0]
HSC_y_lam = y_fil[np.where(y_fil[:,1] == y_fil[:,1].max())[0][0],0]

z_filt_flam = cal_filt_flam(array_spec, z_fil)
plt.scatter(HSC_z_lam, z_filt_flam, c='green')
y_filt_flam = cal_filt_flam(array_spec, y_fil)
plt.scatter(HSC_y_lam, y_filt_flam, c='orange')

obs_z_filt_flam = 10**(-0.4*(m_zAB+2.402+5.0*np.log10(HSC_z_lam)))
obs_y_filt_flam = 10**(-0.4*(m_yAB+2.402+5.0*np.log10(HSC_y_lam)))

plt.scatter(HSC_z_lam, obs_z_filt_flam, c='green', marker="*", s =200, zorder = 10, label='observed z_mag')
plt.scatter(HSC_y_lam, obs_y_filt_flam, c='orange', marker="*", s =200, zorder = 10, label='observed y_mag')

plt.tick_params(labelsize=15)
plt.xlabel("A",fontsize=27)
plt.ylabel(r"f$_\lambda$ erg/s/cm$^2$/A",fontsize=27)
plt.xlim(5000, 20000)
plt.legend(prop={'size':15})

plt.show()

z_mag = -2.5 * np.log10(z_filt_flam ) - 2.402 - 5.0 * np.log10(z_fil[np.where(z_fil[:,1] == z_fil[:,1].max())[0][0],0])
y_mag = -2.5 * np.log10(y_filt_flam ) - 2.402 - 5.0 * np.log10(y_fil[np.where(y_fil[:,1] == y_fil[:,1].max())[0][0],0])
print("z_mag", round(z_mag,3))
print("y_mag", round(y_mag,3))
print("M_1450", round(M1450,3))