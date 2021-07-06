#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 24 16:26:13 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
line_means = ['id','z','mbh','mbh_err','stellar_mass','lbol','spectra','bit','ps_gmag','ps_rmag','ps_imag','ps_rmag','ps_zmag','ps_ymag','host_gmag','host_rmag','host_imag','host_zmag','host_ymag']
infers  = np.loadtxt('HSC_fitting/sdss_quasar_mbh.txt', dtype=str)
IDs_ = infers[:, 0]
HSC_z_ = infers[:,1].astype(np.float)
HSC_Mstar_ = infers[:,4].astype(np.float)
HSC_MBHs_ = infers[:,2].astype(np.float)
HSC_ps_mag = infers[:,10].astype(np.float)
# HSC_ps_mag  =  np.nan_to_num(HSC_ps_mag )
HSC_MBHs_err_ = infers[:,3].astype(np.float)
HSC_label_ = infers[:,-1]
HSC_Lbol = infers[:,5].astype(np.float)


HSC_z = HSC_z_
HSC_Mstar = HSC_Mstar_
HSC_MBHs = HSC_MBHs_

import scipy.stats as st
    
plt.figure(figsize=(11.5,10))      
plt.scatter(HSC_Lbol, HSC_MBHs,c='orange',alpha=0.2,zorder = 1)
from matplotlib import cm 
# contour_cloud(x=HSC_Lbol, y=HSC_MBHs, cmap=cm.Blues)
# plt.scatter(BH_Mass, sdss_g_pointsource)
# plt.plot(np.linspace(5,10), np.linspace(5,10)*0 + abs_Mags)
plt.xlabel('AGN logLbol', fontsize=25)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.tick_params(labelsize=25)
plt.xlim(30, 46.5)
plt.ylim(5.8,10)
xmin, xmax = int(HSC_Lbol.min()-2), int(HSC_Lbol.max()+3)
ymin, ymax = int(HSC_MBHs.min()-2), int(HSC_MBHs.max()+3)
xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([HSC_Lbol, HSC_MBHs])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
plt.contourf(xx, yy, f, cmap=cm.Blues, alpha=0.5)
t = [kernel.pdf([HSC_Lbol[i] , HSC_MBHs[i]]) for i in range(len(HSC_Lbol))]
print(kernel.pdf([45,8]))
