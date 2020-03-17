#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 10:45:35 2020

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
#To estimate the AB mag for a given stellar template
filename  = 'SED_files/summary_*_PA00.fits'
filename = glob.glob(filename)[0]
hdul_sum = pyfits.open(filename)
name_sum = hdul_sum[1].columns
table_sum = hdul_sum[1].data
z_idx = [i for i in range(len(name_sum)) if 'zmc' in str(name_sum[i])][0]
stellar_idx = [i for i in range(len(name_sum)) if 'Msun' in str(name_sum[i])][0]
mel_idx = [i for i in range(len(name_sum)) if 'logZsun' in str(name_sum[i])][0]

z = table_sum[1][z_idx] #Redshift
stellar_ini = table_sum[1][stellar_idx] #stellar mass inferred, to be normed, in M_star
mel = table_sum[1][mel_idx] #Metallicity 

filename_p  = 'SED_files/SFH_*_PA00_param.fits'
filename_p = glob.glob(filename_p)[0]
hdul = pyfits.open(filename_p)
table = hdul[1].data
name = hdul[1].columns
age_idx = [i for i in range(len(name)) if 'T_MW' in str(name[i])][0]
age = str(round(10**table[1][age_idx],3)) # 'AGE:' Gyr
fnu_ini = 0.5 #Input fit always 0.5

#Template used :
print('age:', age, 'mel', round(mel,3), 'sample redshift', z)

def esti_abmag(stellar , fnu_ini = fnu_ini, stellar_ini = stellar_ini):
	'''
	Parameter
	--------
	stellar: In unit of log(M_sun)
	stellar_ini: In unit of M_sun
	'''
	ratio = 10**stellar/stellar_ini
	fnu = fnu_ini * ratio
	zp_fnu = 25
	ab_mag = -2.5*np.log10(fnu) + zp_fnu
	return ab_mag

s_mass = 10.3
print("stellar mass", s_mass, "esti, AB mag:", esti_abmag(s_mass))