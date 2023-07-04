#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 14:19:19 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle, copy, matplotlib
from galight.tools.astro_tools import plt_fits_color
from photutils.aperture import RectangularAperture
from matplotlib.colors import LogNorm   


def scale_bar(ax, d, dist=1/0.13, text='1"', text2=None, color='black', flipped=False, fontsize=20):
    p0 = d / 7.
    ax.plot([p0, p0 + dist], [p0, p0], linewidth=2, color=color)
    ax.text(p0 + dist / 2., p0 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')
    if text2 is not None:
        ax.text(p0 + dist / 2., p0 - 0.08 * d, text2, fontsize=fontsize, color=color, ha='center')
        

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

#!!!
name = ['cid_473', 'cid_1210', 'cid_1245'][0]
image_list = pickle.load(open('colorimage_bin2_{0}.pkl'.format(name),'rb'))
images = []
# zp_list = []
if name != 'cid_1245':
    color_list = [-1,-2,-4]
if name == 'cid_1245':
    color_list = [-1,-3,-4]

filters = [['F115W', 'F150W','F277W', 'F444W'][i] for i in color_list]
use_filt = filters[0]+ '+'+ filters[1]+'+' + filters[2]

for i in color_list:  #['F115W', 'F150W','F277W', 'F444W']
    images.append(image_list[i])
from galight.tools.astro_tools import plt_fits_color, plt_fits
# plt_fits_color(images, Q=7, stretch=0.3)


if name == 'cid_473':
    z = 1.456
if name == 'cid_1210':
    z = 1.924
if name == 'cid_1245':
    z = 0.759

deltaPix = 0.06  #For 150W Bin the 0.015 by factor of 2
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
kpc_per_pixel = scale_relation 
scale = 0.5 * kpc_per_pixel

from astropy.visualization import make_lupton_rgb
rgb_default = make_lupton_rgb(images[0], images[1], images[2], Q=7, stretch=0.3)
fig, ax = plt.subplots()
plt.imshow(rgb_default, origin='lower')
plt.text(1,80,use_filt,fontsize=20, color = 'white')
pixel_s = 0.03*2
ax.set_xticks(np.arange(0,len(images[0]), 1/pixel_s))
ax.set_yticks(np.arange(0,len(images[0]), 1/pixel_s))
ax.set_xticklabels(np.arange(0,len(images[0]), 1/pixel_s)*pixel_s)
ax.set_yticklabels(np.arange(0,len(images[0]), 1/pixel_s)*pixel_s)
scale_bar(ax, len(images[0]), dist=0.5/deltaPix, text='0.5"', text2 ='~{0:.2f}kpc'.format(scale), color = 'white')
plt.savefig('/Users/Dartoon/Downloads/'+name+"_color.pdf",bbox_inches='tight')
plt.show()

#%%
Ebv_image = pickle.load(open('E(BV)_{0}.pkl'.format(name),'rb'))
Av_image = pickle.load(open('Av_{0}.pkl'.format(name),'rb'))
norm = None
# plt.imshow(Ebv_image, norm=norm, origin='lower' ) 
# plt.colorbar()
# plt.show()

f_center = len(Ebv_image)/2

if name == 'cid_473':
    w = 0.2 / deltaPix
    h = 1 / deltaPix
    theta = 20
if name == 'cid_1210':
    w = 0.2 / deltaPix
    h = 1 / deltaPix
    theta = -20
if name == 'cid_1245':
    w = 0.2 / deltaPix
    h = 1.5 / deltaPix
    theta = -20

aper = RectangularAperture((f_center, f_center), w, h, theta=theta/180*np.pi)
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat'))  # copy the default cmap
my_cmap.set_bad('black')
vmin = 1.e-3
# plt.imshow(Ebv_image, origin='lower', cmap=my_cmap, norm=LogNorm())  # , vmin=vmin, vmax=vmax)
plt.imshow(Ebv_image, norm=norm, origin='lower',vmin = vmin) 
plt.colorbar()
aper.plot(color='blue',
          lw=3.5)
plt.show()
from photutils.aperture import aperture_photometry
mask_ = aper.to_mask(method = 'center')
mask = mask_.to_image(Ebv_image.shape)
mean = np.mean(Ebv_image[mask==1])
std = np.std(Ebv_image[mask==1])
print('{2} E(B-V): {0:.3}$\pm${1:.3}'.format(mean,std, name))
print("NO. of pixel ", np.sum(mask==1))

mean = np.mean(Av_image[mask==1])
std = np.std(Av_image[mask==1])
print('{2} Av: {0:.3}$\pm${1:.3}'.format(mean,std, name))
print("NO. of pixel ", np.sum(mask==1))
# phot_table_host = to_mask(Ebv_image, aper)

mean = np.mean(Av_image[mask==1]/4.05)
std = np.std(Av_image[mask==1]/4.05)
print('{2} Av/4.05: {0:.3}$\pm${1:.3}'.format(mean,std, name))
print("NO. of pixel ", np.sum(mask==1))
# phot_table_host = to_mask(Ebv_image, aper)

#%%Filter flux
zp = 28.0866
filt_list = ['F115W', 'F150W','F277W', 'F444W']
print(name,'zeropoint:',zp)
# print('zeropoint:',zp)
for i, filt in enumerate(filt_list):
    image = image_list[i]
    mean = np.mean(image[mask==1])
    std = np.std(image[mask==1])
    print(filt,round(mean,4),'$\pm$',round(std,4))
# 
#%%
EBV = Av_image/3.1
factor = 5
bools = (images[2] < np.std(images[2][:,:3])*factor) * (images[1] < np.std(images[1][:,:3])*factor ) *  (images[0] < np.std(images[0][:,:3])*factor) 
EBV[ bools ] = 0 
print(name, "overall mean E(B-V):", np.mean(EBV[ EBV !=0 ] ) )
fig, ax = plt.subplots()
plt.imshow(EBV, norm=norm, origin='lower',vmin = vmin) 
# aper.plot(color='red',
#           lw=1.5,alpha = 0.7)
plt.colorbar()
pixel_s = 0.03*2
ax.set_xticks(np.arange(0,len(images[0]), 1/pixel_s))
ax.set_yticks(np.arange(0,len(images[0]), 1/pixel_s))
ax.set_xticklabels(np.arange(0,len(images[0]), 1/pixel_s)*pixel_s)
ax.set_yticklabels(np.arange(0,len(images[0]), 1/pixel_s)*pixel_s)
scale_bar(ax, len(images[0]), dist=0.5/deltaPix, text='0.5"', text2 ='{0:.2f}kpc'.format(scale), color = 'white')
plt.savefig('/Users/Dartoon/Downloads/'+name+"_E(B-V).pdf")
plt.show()




# Ebv_image_scale = Ebv_image/Ebv_image.max()*Av_image.max()/4.05
# Ebv_image_scale[bools] = 0
# plt.imshow(Ebv_image_scale, norm=norm, origin='lower',vmin = vmin) 
# plt.tick_params(labelsize=15)
# plt.colorbar()
# plt.savefig('/Users/Dartoon/Downloads/'+name+"E(B-V).pdf")
# plt.show()
