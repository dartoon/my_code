#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 20:48:49 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'


#shibuya2015, lg(r0/kpc) =-0.38+ -0.56,  r0 = 0.419+-1.12  beta = 0.25+-0.05

plt.figure(figsize=(11, 7))

#Target 0
Muv_0 = [-18.35, -18.46, -18.26]
R_filt_0 = np.array([1.87, 1.06 ])
q = 0.630
R_filt_0 = R_filt_0*np.sqrt(q)
# wave = F356W
smass = 10.63
z = 6.34
key = 'F356W'
filt_wave_dic = {'F115W': 11623.89, 'F150W': 15104.23, 'F200W': 20028.15, 
           'F277W': 27844.64, 'F356W': 35934.49, 'F444W': 44393.50, 'F410M': 40886.55}
filt_wave = filt_wave_dic[key]
wave_rest = 1600
wave = filt_wave / (1+z)
z_p = int(key[1:4])/(wave_rest/100)  -1
Reff_filt = R_filt_0[0]
vdep = -0.35 + 0.12 * z - 0.25 * np.log10( 10 ** smass / 10**10)
R_0 = [Reff_filt * ((1+z)/(1+z_p)) ** vdep]
R_0.append(R_0[0]* R_filt_0[1]/R_filt_0[0] )
plt.scatter(Muv_0[0], np.log10(R_0[0]),c='blue',s=280,marker="o",
            zorder=0, edgecolors='black',alpha=1, label = 'J2255+0251')
plt.errorbar(Muv_0[0], np.log10(R_0[0]), 
             yerr = [ [np.log10(R_0[0]) - np.log10(R_0[0]-R_0[1] )], 
                     [ np.log10(R_0[0]+R_0[1] ) - np.log10(R_0[0])] ],
              # xerr = [ [Muv_0[0] - Muv_0[1]], [Muv_0[2] - Muv_0[0]] ],
              xerr = [ [0.6], [0.6] ],
              c='blue', elinewidth=3)

#Target 1
Muv_0 = [-20.57, -20.74, -20.4]
key = 'F150W'
# key = 'F356W'
if key == 'F356W':
    R_filt_0 = np.array([0.75, 0.11])
    q = 0.360
if key == 'F150W':
    R_filt_0 = np.array([0.59, 0.14])
    q = 0.269
R_filt_0 = R_filt_0*np.sqrt(q)
smass = 11.14
z = 6.4
filt_wave = filt_wave_dic[key]
wave = filt_wave / (1+z)
z_p = int(key[1:4])/(wave_rest/100)  -1
Reff_filt = R_filt_0[0]
vdep = -0.35 + 0.12 * z - 0.25 * np.log10( 10 ** smass / 10**10)
R_0 = [Reff_filt * ((1+z)/(1+z_p)) ** vdep]
R_0.append(R_0[0]* R_filt_0[1]/R_filt_0[0] )
plt.scatter(Muv_0[0], np.log10(R_0[0]),c='orange',s=280,marker="o",
            zorder=10, edgecolors='black',alpha=1, label = 'J2236+0032')
plt.errorbar(Muv_0[0], np.log10(R_0[0]), 
             yerr = [ [np.log10(R_0[0]) - np.log10(R_0[0]-R_0[1] )], 
                     [ np.log10(R_0[0]+R_0[1] ) - np.log10(R_0[0])] ],
              # xerr = [ [Muv_0[0] - Muv_0[1]], [Muv_0[2] - Muv_0[0]] ],
              xerr = [ [0.6], [0.6] ], c='orange', zorder=10,
              elinewidth=3)

#%%
M=np.linspace(-24,-17,10)
r0_s = [-0.2, 0.25]
beta_s = [0.27, 0.01]
r0_liter = [r0_s]
beta_liter = [beta_s]
color_liter = ['black']
label_liter = [' Shibuya+2015']
k=0
r0sim=r0_liter[k][0] + r0_liter[k][1] * np.random.randn(100)
bsim=beta_liter[k][0] + beta_liter[k][1] * np.random.randn(100)
X=[]
for i in range(100):
    X.append(-bsim[i]*0.4*(M+21)+r0sim[i])
X=np.array(X)
mu=X.mean(axis=0)
sigma=X.std(axis=0)
plt.plot(M, mu, lw=5, color=color_liter[k], label =label_liter[k],zorder=3)
plt.fill_between(M, (mu+sigma), (mu-sigma), facecolor=color_liter[k], alpha=0.15)
#%%Plot Shibuya:
f = open("../other_files/shibuya","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
lines = [line for line in lines if 'z6' in line or 'z7' in line]

import astropy.units as u
from astropy.cosmology import LambdaCDM, FlatLambdaCDM
cosmo1 = LambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, 0.7)

arc_per_kpc = {6:cosmo1.arcsec_per_kpc_proper(6).value, 7: cosmo1.arcsec_per_kpc_proper(7).value}
distmod = {6:cosmo1.distmod(6).value, 7: cosmo1.distmod(7).value}
# 5*np.log10(cosmo1.luminosity_distance(z=6).value * 10**6) -5

shibuya_scatter = []
for line in lines:
    z = float(line[1])
    info = line.split(' ')
    info = [float(info[i]) for i in range(1,len(info)) if '' != info[i]]
    muv = info[0]
    emuv = info[1]
    Muv = muv - distmod[z] + 2.5 * np.log10(1+z) 
    reff = info[2] * np.sqrt(info[4]) / arc_per_kpc[z] # info[4] is q.
    ereff = info[3] * np.sqrt(info[4])  / arc_per_kpc[z] # info[4] is q.
    shibuya_scatter.append( [Muv, emuv, reff, ereff])
shibuya_scatter = np.array(shibuya_scatter)
plt.scatter(shibuya_scatter[:,0], np.log10(shibuya_scatter[:,2]),c='gray',s=180,marker=".",
            zorder=0, edgecolors='black',alpha=1, label = 'Shibuya')
# plt.errorbar(Muv_0[0], np.log10(R_0[0]), 
#              yerr = [ [np.log10(R_0[0]) - np.log10(R_0[0]-R_0[1] )], 
#                      [ np.log10(R_0[0]+R_0[1] ) - np.log10(R_0[0])] ],
#               # xerr = [ [Muv_0[0] - Muv_0[1]], [Muv_0[2] - Muv_0[0]] ],
#               xerr = [ [0.6], [0.6] ],
#               c='blue', elinewidth=3)

#%%
plt.xlabel(r"M$_{\rm UV}$",fontsize =30)#,weight='bold')
plt.ylabel(u'log($R_{\mathrm{eff, circ.}}$, units of kpc)',fontsize =30)#,weight='bold')
plt.tick_params(labelsize=20)
plt.ylim([-1.5, 1.5])
plt.xlim([-23, -17])
plt.legend(scatterpoints=1,numpoints=1,loc=1,prop={'size':20},ncol=2,handletextpad=0)
plt.savefig('figures/z6-7_size_luminosity.pdf')
plt.show()