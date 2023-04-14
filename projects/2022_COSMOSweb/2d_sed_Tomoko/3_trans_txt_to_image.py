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

def scale_bar(ax, d, dist=1/0.13, text='1"', text2=None, color='black', flipped=False, fontsize=20):
    p0 = d / 7.
    ax.plot([p0, p0 + dist], [p0, p0], linewidth=2, color=color)
    ax.text(p0 + dist / 2., p0 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')
    if text2 is not None:
        ax.text(p0 + dist / 2., p0 - 0.08 * d, text2, fontsize=fontsize, color=color, ha='center')
        

f = open('fmos_alma_cosmosweb.cat','r')
string = f.read()
lines = string.split('\n')
lines = [lines[i] for i in range(len(lines)) if 'FMOS_J09' in lines[i]]
idx = 4

target_name, RA, Dec, z, best_mass = lines[idx].split(' ')
name = target_name[7:12]
t_name = target_name[5:12]
 
sed_2d_info = pickle.load(open('2d_filts_mag_bin2_{0}.pkl'.format(t_name),'rb'))

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
        count, smass, sfr, m_age, l_age, AV = line.split(' ')
        count = int(count[4:])
        _i, _j = sed_2d_info[count][0], sed_2d_info[count][1]
        smass_image[_i, _j] = float(smass)    #smass in logMsun
        sfr_image[_i, _j] = float(sfr)          #logMsun/yr 
        age_image[_i, _j] = 10**float(m_age)    #logGyr to Gry
        AV_image[_i, _j] = AV    #logGyr
        # Ebv_image[_i, _j] = Ebv    #E(B-V)
        # Rv.append(float(AV)/float(Ebv))

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
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
from galight.tools.astro_tools import plt_fits_color
from astropy.cosmology import FlatLambdaCDM
from astropy.visualization import make_lupton_rgb
image_list = pickle.load(open('colorimage_bin2_{0}.pkl'.format(t_name),'rb'))

images = []
use_filt = ''
use_filt_id = [4, 3, 1]
filts = ['F814W','F115W', 'F150W','F277W', 'F444W']
for i in use_filt_id:  #['ACS','F115W', 'F150W','F277W', 'F444W']
    images.append(image_list[i])
    use_filt += filts[i] + '+'
use_filt = use_filt[:-1]
deltaPix = 0.03*2
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
kpc_per_pixel = scale_relation 
scale = 0.5 * kpc_per_pixel
rgb_default = make_lupton_rgb(images[0], images[1], images[2], Q=7, stretch=0.3)
fig, ax = plt.subplots()
plt.title(t_name)
plt.imshow(rgb_default, origin='lower')
plt.text(1,80,use_filt,fontsize=20, color = 'white')
scale_bar(ax, len(images[0]), dist=0.5/deltaPix, text='0.5"', text2 ='~{0:.2f}kpc'.format(scale), color = 'white')
plt.savefig('products/color_image_{0}.pdf'.format(t_name))
plt.show()


norm = None  
print('smass_image')
norm = LogNorm(vmin=4.5, vmax=8)#np.max(img[~np.isnan(img)]))
plt.title(t_name)
plt.imshow(smass_image, norm=norm, origin='lower' ) 
plt.colorbar()
plt.show()
print(np.log10(np.sum(10**smass_image)))
print("BestMass", best_mass)

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
plt.title(t_name)
plt.imshow(AV_image, norm=norm, origin='lower' ) 
plt.colorbar()
plt.show()
pickle.dump(AV_image , open('Av_{0}.pkl'.format(name), 'wb'))

pyfits.PrimaryHDU(smass_image).writeto('products/2D_logsmass_{0}.fits'.format(t_name),overwrite=True)
pyfits.PrimaryHDU(AV_image).writeto('products/2D_Av_{0}.fits'.format(t_name),overwrite=True)

# import matplotlib
# cmap_r = matplotlib.cm.get_cmap('RdBu_r')
