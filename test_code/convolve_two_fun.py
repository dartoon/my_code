#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 16:07:43 2019

@author: dxh

test convolving two func
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from scipy.integrate import quad

import matplotlib as matt
import matplotlib.lines as mlines
from matplotlib import colors
matt.rcParams['font.family'] = 'STIXGeneral'

def poss_ln_gaussian(m1, mu, sigma):
    if m1 >= 0:
        poss =  1/(m1 * sigma * np.sqrt(2 * np.pi)) * np.exp(-(np.log(m1) - mu)**2/(2*sigma**2))
    else:
        poss = 0
    return poss


def mass_fun_i(m1, a=2.35, mbh_max=80, mbh_min=5):
    if m1>=mbh_min and m1<=mbh_max:
        N = (mbh_max)**(-a+1)/(1-a) - (mbh_min)**(-a+1)/(1-a)             # The intergral function of m**a for normalize the power law function
        return m1**(-a)/N
    else:
        return 0.
    
poss_mass_fun=np.vectorize(mass_fun_i)    

def joint_twofun(tau, t, a=2.35, mbh_max=80, mbh_min=5, sigma =0.2):
    return mass_fun_i(tau, a=a, mbh_max=mbh_max, mbh_min=mbh_min) * poss_ln_gaussian(t, mu=np.log(tau), sigma = sigma)
def cov_twofun(t, a=2.35, mbh_max=80, mbh_min=5, sigma =0.2):
    inter = quad(joint_twofun, 0, 100, args=(t, a, mbh_max, mbh_min,sigma))[0]
    return inter
cov_twofun=np.vectorize(cov_twofun)   


plt.subplots(figsize=(8, 6))
#print cov_twofun(5.5)
#Before conv:
x = np.logspace(np.log10(0.2), np.log10(85),300)
y = poss_mass_fun(x)
#y = poss_ln_gaussian(x,mu=np.log(5.5),sigma=0.2)
plt.plot(x, y ,'--', label = 'power-law distribution')
#plt.xlim(3,20)
#plt.show()

#After conv:
x = np.logspace(np.log10(0.2), np.log10(85),300)
#y = poss_mass_fun(x)
#y = poss_ln_gaussian(x,mu=np.log(5.5),sigma=0.2)
y = cov_twofun(x, sigma=0.04)
plt.plot(x, y, label = 'power-law convolved with Log-Normal')
plt.xlabel("$m_1 (M_{\odot})$",fontsize=25)
plt.ylabel("P($m$)",fontsize=25)
plt.xlim(3,20)
plt.tick_params(labelsize=15)
plt.legend(prop={'size': 17},loc=1)
plt.savefig("convolving.pdf")
plt.show()
