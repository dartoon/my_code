#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 10:20:29 2021

@author: Dartoon
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

# Coumn 1: redshift
# Coumn 2: log stellar mass
# Coumn 3: log BH mass
# Coumn 4: g-band magnitude of host galaxy 
# Coumn 5: r-band magnitude of host galaxy
# Coumn 6: AGN bolometric luminosity in units of 10^45 erg/s
# Coumn 7: statistical weight
filename = 'SAM/catalogue.dat'
text  = np.loadtxt(filename)
text = text[text[:,5] != 0]
# print(text.shape)
redshift = text[:, 0]
logM_mass =  text[:, 1]
logMBH =  text[:, 2]
g_galaxy =  text[:, 3]
r_galaxy =  text[:, 4]
AGN_bol_10_45 =  text[:, 5]
stat_weight =  text[:, 6]
# plt.scatter(logM_mass, g_galaxy)  #Look at the scatter.

L_5100_10_45 = AGN_bol_10_45 / 9.26
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dl = cosmo.luminosity_distance(redshift).value  # Mpc
dis = dl * 3.085677581e+24 #from Mpc to cm
F_lam_5100 = L_5100_10_45*10**45 / (4*np.pi * dis**2) / 5100
wave_lam_g = 4700 * (1+redshift)  #rest-frame in A
F_lam_g = F_lam_5100 / 5100**(-(-0.44+2)) * wave_lam_g **(-(-0.44+2)) #erg/s/cm^2/A
F_mu_g = F_lam_g *  (wave_lam_g) **2 / (2.9979 * 10**18)
obs_mag_g = [(-2.5 * np.log10(F_mu_g[i]) - 48.60) for i in range(len(F_mu_g))]
abs_mag_g = obs_mag_g - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
wave_lam_r = 6100 * (1+redshift)  #rest-frame in A  
F_lam_r = F_lam_5100 / 5100**(-(-0.44+2)) * wave_lam_r **(-(-0.44+2)) #erg/s/cm^2/A
F_mu_r = F_lam_r *  (wave_lam_r) **2 / (2.9979 * 10**18)
obs_mag_r = [(-2.5 * np.log10(F_mu_r[i]) - 48.60) for i in range(len(F_mu_r))]
abs_mag_r = obs_mag_r - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
