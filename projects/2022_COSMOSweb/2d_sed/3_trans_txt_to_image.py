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


name = ['cid_473', 'cid_1210', 'cid_1245'][2]
 
sed_2d_info = pickle.load(open('2d_filts_mag_bin2_{0}.pkl'.format(name),'rb'))

f = open("{0}_sed_2d_result_bin2.txt".format(name),"r")
string = f.read()
lines = string.split('\n')   # Split in to \n
size = int(np.sqrt(len(sed_2d_info)))
smass_image = np.zeros([size,size])
sfr_image =  np.zeros([size,size])
age_image =  np.zeros([size,size])
AV_image = np.zeros([size,size])
Ebv_image = np.zeros([size,size])

Rv = []
for ct, line in enumerate(lines[1:-1]):
    if len(line.split(' ')) < 4:
        continue
    else:
        count, smass, sfr, m_age, l_age, AV, Ebv = line.split(' ')
        count = int(count)
        _i, _j = sed_2d_info[count][0], sed_2d_info[count][1]
        smass_image[_i, _j] = float(smass)    #smass in logMsun
        sfr_image[_i, _j] = float(sfr)          #logMsun/yr 
        age_image[_i, _j] = 10**float(m_age)    #logGyr to Gry
        AV_image[_i, _j] = AV    #logGyr
        Ebv_image[_i, _j] = Ebv    #E(B-V)
        Rv.append(float(AV)/float(Ebv))

# for i in range(len(smass_image)):
#     for j in range(len(smass_image)):
#         if smass_image[i,j]>8.:
#             check = np.average([ smass_image[i-1,j], smass_image[i+1,j], 
#                                            smass_image[i,j-1], smass_image[i,j+1]])
#             if smass_image[i,j] > check*1.1:
#                 smass_image[i,j] = check
#                 sfr_image[i,j] = np.average([ sfr_image[i-1,j], sfr_image[i+1,j], 
#                                                sfr_image[i,j-1], sfr_image[i,j+1]])

        
#%%
import pickle 
from galight.tools.astro_tools import plt_fits_color
image_list = pickle.load(open('colorimage_bin2_{0}.pkl'.format(name),'rb'))
images = []
# zp_list = []
for i in [-1,-3,-4]:  #['F115W', 'F150W','F277W', 'F444W']
    # zp_list.append(zp_dict[filts[i]])
    images.append(image_list[i])
from galight.tools.astro_tools import plt_fits_color, plt_fits
# images = [images[i] * 10 ** (-0.4*(zp_list[i]-zp_list[0])) for i in range(len(images)) ]
# for i in range(len(images)):
#     plt_fits(images[i], vmin=0.001, vmax=2.5)
plt_fits_color(images, Q=7, stretch=0.3)


norm = None  
print('smass_image')
norm = LogNorm(vmin=4.5, vmax=8)#np.max(img[~np.isnan(img)]))
plt.imshow(smass_image, norm=norm, origin='lower' ) 
plt.colorbar()
plt.show()

# print('sfr_image')
# norm = LogNorm(vmin=0.003, vmax=0.1)#np.max(img[~np.isnan(img)]))
# plt.imshow(sfr_image, origin='lower' ) 
# plt.colorbar()
# plt.show()


# print('age_image')
# norm = LogNorm(vmin=0.002, vmax=3)#np.max(img[~np.isnan(img)]))
# plt.imshow(age_image, norm=norm, origin='lower' ) 
# plt.colorbar()
# plt.show()


# print('AV_image')
# # norm = LogNorm(vmin=0.001, vmax=3)#np.max(img[~np.isnan(img)]))
# norm = None
# plt.imshow(AV_image, norm=norm, origin='lower' ) 
# plt.colorbar()
# plt.show()

# print('Ebv_image')
# # norm = LogNorm(vmin=0.001, vmax=3)#np.max(img[~np.isnan(img)]))
# norm = None
# plt.imshow(Ebv_image, norm=norm, origin='lower' ) 
# plt.colorbar()
# plt.show()
# pickle.dump(Ebv_image , open('E(BV)_{0}.pkl'.format(name), 'wb'))

print('Av_image')
# norm = LogNorm(vmin=0.001, vmax=3)#np.max(img[~np.isnan(img)]))
norm = None
plt.imshow(AV_image, norm=norm, origin='lower' ) 
plt.colorbar()
plt.show()
pickle.dump(AV_image , open('Av_{0}.pkl'.format(name), 'wb'))



# import matplotlib
# cmap_r = matplotlib.cm.get_cmap('RdBu_r')
