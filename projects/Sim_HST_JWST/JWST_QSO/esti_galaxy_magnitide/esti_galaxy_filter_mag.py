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

ID = 1 #1, 2, 3, 4, 5
#The data information of the five system
if ID ==1:
    M1450, m_zAB, m_yAB, z, M_dyn = -23.82, 22.775, 22.942, 6.10, 82.
elif ID ==2:    
    M1450, m_zAB, m_yAB, z, M_dyn = -25.31, 21.827, 21.614, 6.37, 14.
elif ID ==3:
    M1450, m_zAB, m_yAB, z, M_dyn = -24.09, 22.768, 23.649, 6.39, 56.
elif ID ==4:    
    M1450, m_zAB, m_yAB, z, M_dyn = -24.73, 22.121, 22.045, 6.2, 13.
elif ID ==5:
    M1450, m_zAB, m_yAB, z, M_dyn = -24.69, 22.402, 22.326, 6.26, 290.

def esti_abmag(stellar , fnu_ini, stellar_ini):
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

#%%
filt_l = ['444','356', '200', '150']
#To estimate the AB mag for a given stellar template
for i in range(len(filt_l)):
    filt = filt_l[i]
    filename  = 'f{0}w/summary_{0}{1}_PA00.fits'.format(filt,ID)
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
    filename_p  = 'f{0}w/SFH_{0}{1}_PA00_param.fits'.format(filt,ID)
    filename_p = glob.glob(filename_p)[0]
    hdul = pyfits.open(filename_p)
    table = hdul[1].data
    name = hdul[1].columns
    age_idx = [i for i in range(len(name)) if 'T_MW' in str(name[i])][0]
    age = str(round(10**table[1][age_idx],3)) # 'AGE:' Gyr
    fnu_ini = 0.5 #Input fit always 0.5
    #Template used :
    if i ==0:
        print("ID", ID, 'age:', age, 'mel', round(mel,3), 'sample redshift', z) #Same for all the filters
    #s_mass = 10.3
    s_mass = np.log10(M_dyn * 10 **9) 
    print("\nintput stellar mass:", round(s_mass,3), "\nesti, F{0}W AB mag:".format(filt), round(esti_abmag(s_mass, fnu_ini, stellar_ini),3))
