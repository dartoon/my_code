#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 20:33:39 2021

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob

filenames = glob.glob('TNG100/*') 
filenames.sort()

# for i in range(len(filenames)):
for i in [4]:    
    text  = np.load(filenames[i])
    zs = float(filenames[i].split("_z")[1][:4])
    print(text.shape)
Stellar_Mass, BH_Mass, sdss_i_galaxy, sdss_g_galaxy, sdss_r_galaxy, sdss_i_pointsource, sdss_g_pointsource = np.load('TNG100/data_z0.40.npy')
plt.scatter(np.log10(Stellar_Mass), sdss_g_galaxy)

#%%
I_mag_break = 21.7

import numpy as np
h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
from scipy.integrate import quad
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)               #Perform the function EE in array style

#Calculate the luminosity distance:
dl=(1+zs)*c*vec_EE(zs)/h0 *10**6     #zs is a list of redshift of the sources (input as 1 D numpy array)

#Transfer the magnitude to absolute ones:
abs_Mags = I_mag_break -5*(np.log10(dl)-1)   #dl is the luminosity distance which is a function of redshift:

    
sdss_g_totall = -2.5*np.log10(10**(-0.4*sdss_g_galaxy) + 10**(-0.4*sdss_g_pointsource))

plt.figure(figsize=(11.5,12))
Stellar_Mass = Stellar_Mass[sdss_g_totall<abs_Mags]
BH_Mass = BH_Mass[sdss_g_totall<abs_Mags]

BH_Mass = np.log10(BH_Mass)
Stellar_Mass = np.log10(Stellar_Mass)

dMBH, dmag, dMstar= 0.4, 0.3, 0.17  #dmag is for host magnitude. 

Stellar_Mass_nois = Stellar_Mass + np.random.normal(0, dMstar, size=Stellar_Mass.shape)
BH_Mass_nois = BH_Mass + np.random.normal(0, dMBH, size=BH_Mass.shape)


plt.scatter(Stellar_Mass_nois, BH_Mass_nois,c='red',
            s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.4)

plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
plt.xlim(9,12.5)
plt.ylim(6.0,10.3)
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)
plt.show()
