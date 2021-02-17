#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 16:38:02 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
infers  = np.loadtxt('./local_Vardha_data.txt', dtype=str)
ID = infers[:, 0]
z = infers[:, 3].astype(np.float)
MBH = infers[:, 6].astype(np.float)
Mstar = infers[:, 11]  #sph mass
VD = infers[:, 7].astype(np.float)

for i in range(len(Mstar)):
    if Mstar[i] == '...':
        Mstar[i] = -99
    else:
        Mstar[i] = float(Mstar[i].split('$\\pm$')[0])
Mstar = Mstar.astype(np.float)

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
da = cosmo.angular_diameter_distance(z).value  # Mpc

infers_t3  = np.loadtxt('./local_Vardha_data_t3.txt', dtype=str)
ID_t3 = infers_t3[:, 0]
Re_sph = infers_t3[:, -3].astype(np.float)
scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
Re_sph_kpc = Re_sph * scale_relation   #In arcsec

#%%Load Kormendy 2013:
kormendy_t2 = np.loadtxt('./data/kormendy13_t2.txt', dtype=str)

kor_dis = kormendy_t2[:, 0].astype(np.float)
kor_Mstar = kormendy_t2[:, 6]
kor_sig = kormendy_t2[:, 9]
kor_Reff = kormendy_t2[:, ] #Could find such information...

#%%

#Plot the fundamental plan elements:
logR = np.log10(Re_sph_kpc)
a = 1.629
b = -0.84
right = a * np.log10(VD) + b * np.log10(10**Mstar/(2*np.pi*Re_sph_kpc**2) )
gam = 4.496

plt.figure(figsize=(10, 10))
plt.plot(np.linspace(-0.8, 0.6), np.linspace(-0.8, 0.6) - gam ,'black', label = 'Bezanson2015, total')
plt.scatter(logR, right, label = 'Vardha 2021')
plt.ylim(-5.8, -3.3)
# plt.xlim(-0.5,1.5)
plt.xlabel("log Reff [kpc]",fontsize=27)
plt.ylabel("1.63 log $\sigma$ - 0.84 log$\Sigma_*$ ",fontsize=27)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':28})
plt.show()