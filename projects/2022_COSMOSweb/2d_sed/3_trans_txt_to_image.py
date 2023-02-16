#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 11:19:46 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

#     # host_residual_list = fit_run.
import pickle
from matplotlib.colors import LogNorm     


name = ['cid_473', 'cid_1210', 'cid_1245'][0]
 
sed_2d_info = pickle.load(open('2d_filts_mag_bin2_{0}.pkl'.format(name),'rb'))
f = open("{0}_sed_2d_result_bin2.txt".format(name),"r")
string = f.read()
lines = string.split('\n')   # Split in to \n
size = int(np.sqrt(len(sed_2d_info)))
smass_image = np.zeros([size,size])
sfr_image =  np.zeros([size,size])
age_image =  np.zeros([size,size])
AV_image = np.zeros([size,size])
for ct, line in enumerate(lines[1:-1]):
    if len(line.split(' ')) < 4:
        continue
    else:
        count, smass, sfr, m_age, l_age, AV = line.split(' ')
        count = int(count)
        _i, _j = sed_2d_info[count][0], sed_2d_info[count][1]
        smass_image[_i, _j] = float(smass)    #smass in logMsun
        sfr_image[_i, _j] = float(sfr)          #logMsun/yr 
        age_image[_i, _j] = 10**float(m_age)    #logGyr to Gry
        AV_image[_i, _j] = AV    #logGyr

        
for i in range(len(smass_image)):
    for j in range(len(smass_image)):
        if smass_image[i,j]>8.:
            check = np.average([ smass_image[i-1,j], smass_image[i+1,j], 
                                           smass_image[i,j-1], smass_image[i,j+1]])
            if smass_image[i,j] > check*1.1:
                smass_image[i,j] = check
                sfr_image[i,j] = np.average([ sfr_image[i-1,j], sfr_image[i+1,j], 
                                               sfr_image[i,j-1], sfr_image[i,j+1]])

        
#%%

norm = None  
norm = LogNorm(vmin=4.5, vmax=8)#np.max(img[~np.isnan(img)]))
plt.imshow(smass_image, norm=norm, origin='lower' ) 
plt.colorbar()
plt.show()

norm = LogNorm(vmin=0.003, vmax=0.1)#np.max(img[~np.isnan(img)]))
plt.imshow(sfr_image, origin='lower' ) 
plt.colorbar()
plt.show()


norm = LogNorm(vmin=0.002, vmax=3)#np.max(img[~np.isnan(img)]))
plt.imshow(age_image, norm=norm, origin='lower' ) 
plt.colorbar()
plt.show()


# norm = LogNorm(vmin=0.001, vmax=3)#np.max(img[~np.isnan(img)]))
norm = None
plt.imshow(AV_image, norm=norm, origin='lower' ) 
plt.colorbar()
plt.show()


import matplotlib
cmap_r = matplotlib.cm.get_cmap('RdBu_r')
