#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 17:20:27 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from galight.tools.astro_tools import plt_fits

fitsFile = pyfits.open('J2255p0251/jw01967-o012_s00002_nirspec_f290lp-g395m-s200a2_s2d.fits')

spec_2d= fitsFile[1].data
header = fitsFile[0].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']

rx, ry = 27, 33
plt_fits( spec_2d[rx:ry ,400:500] )
plt.show()

#%%
cut1, cut2 = 444,464 #OIII
# cut1, cut2 = 1020,1150 #Halhpa
# cut1, cut2 = 600,900 #Contumim
# cut1, cut2 = 744,964
cut1_r, cut2_r = 474,494 #Contumim
cut1_l, cut2_l = 420,440 #Contumim
x_range = len(spec_2d.T)
plt.plot( np.linspace(0,x_range-1,x_range), np.sum(spec_2d[rx:ry ,:], axis=0))
plt.axvline(x=cut1, ls='--', linewidth=1.6, c='red', zorder = 1, label='true value')
plt.axvline(x=cut2, ls='--', linewidth=1.6, c='red', zorder = 1, label='true value')
plt.axvline(x=cut1_l, ls='--', linewidth=1.6, c='green', zorder = 1, label='true value')
plt.axvline(x=cut2_l, ls='--', linewidth=1.6, c='green', zorder = 1, label='true value')
plt.axvline(x=cut1_r, ls='--', linewidth=1.6, c='green', zorder = 1, label='true value')
plt.axvline(x=cut2_r, ls='--', linewidth=1.6, c='green', zorder = 1, label='true value')
plt.xlim([100,800])
plt.show()

OIII_slit_spec = np.sum(spec_2d[rx:ry,cut1:cut2], axis=1) 
Cont_slit_spec = (np.sum(spec_2d[rx:ry,cut1_l:cut2_l], axis=1)  + np.sum(spec_2d[rx:ry,cut1_r:cut2_r], axis=1) )/2
x_range = len(OIII_slit_spec)
fluxes = OIII_slit_spec

plt.plot( np.linspace(0,(x_range-1)*0.1,x_range), fluxes, 'red', label='OIII region, entire')
plt.plot( np.linspace(0,(x_range-1)*0.1,x_range), Cont_slit_spec, 'green', label='green') #/np.sum(Cont_slit_spec)*np.sum(fluxes))

pure_OIII = fluxes - Cont_slit_spec
plt.plot( np.linspace(0,(x_range-1)*0.1,x_range), pure_OIII, label='Pure OIII line') #/np.sum(Cont_slit_spec)*np.sum(fluxes))
plt.xlabel('arcsec')
plt.legend()
plt.show()


from scipy.optimize import curve_fit
x_data = np.linspace(0,(x_range-1)*0.1,x_range)

def func(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

parameters, covariance = curve_fit(func, x_data, fluxes - Cont_slit_spec)
sigma = abs(parameters[-1])
print("pureOIII:",sigma)
parameters, covariance = curve_fit(func, x_data, Cont_slit_spec)
sigma = abs(parameters[-1])
print("Borad:",sigma)

# #%%
# x_range = len(pure_OIII)
# plt.plot( -0.1+np.linspace(0,(x_range-1)*0.1,x_range), pure_OIII/np.max(pure_OIII), 'red', label='OIII region, NIRSpec')
# x_range = len(read)
# plt.plot( np.linspace(0,(x_range-1)*0.1,x_range), read/read.max(), 'blue', label='NIRCam')
# plt.legend()
# plt.show()
