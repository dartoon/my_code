#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 13:23:13 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from functions_for_result import esti_smass, load_prop, load_info
from scipy.ndimage import zoom

# ID, mags, z = 'idx0', 
# 1,2,0,51,35
idx = 1
# root_folder = '../*/*'  #Include HST
root_folder = './*'  #Not include HST
fit_run_dict = load_prop(idx, root_folder = root_folder, prop_name='fit_run')
filt_list = list(fit_run_dict.keys())
# host_residual_list = []
deltaPix_list = []
#%%
shift_center = True
l_band = 'F356W'
fit_run = fit_run_dict[l_band]
l_deltaPix = fit_run.fitting_specify_class.deltaPix
# pos = fit_run.final_result_ps[0]['ra_image'], fit_run.final_result_ps[0]['dec_image'] 

from galight.tools.astro_tools import plt_fits

run_filt_list = [l_band] + [filt_list[i] for i in range(len(filt_list)) if filt_list[i] != l_band]


image_list = [None] * len(run_filt_list)
for i, filt in enumerate(run_filt_list):
    fit_run = fit_run_dict[filt]
    img_org = fit_run.flux_2d_out['data-Point Source'] - np.sum(fit_run.image_host_list[1:],axis=0 )  #!!!
    deltaPix = fit_run.fitting_specify_class.deltaPix
    # img_org = fit_run.flux_2d_out['model']
    if shift_center == True:
        if hasattr(fit_run, 'final_result_ps'):
            pos =  - fit_run.final_result_ps[0]['ra_image'][0]/deltaPix ,\
                  fit_run.final_result_ps[0]['dec_image'][0]/deltaPix
        else:
            pos =  - fit_run.final_result_galaxy[0]['center_x']/deltaPix ,\
                  fit_run.final_result_galaxy[0]['center_y']/deltaPix
        pos = np.int0(np.array(pos))
        if filt == l_band:
            l_pos = pos
        if filt != l_band:
            shift = pos - l_pos
            new_pos = pos + shift + int(len(img_org)/2) #!!!
            ct = int(len(img_org)/2) - np.max(shift) - 1 
            img_org =  img_org[ new_pos[1] - ct:new_pos[1] + ct+1, new_pos[0] - ct:new_pos[0] + ct +1 ]
    ratio = fit_run.fitting_specify_class.deltaPix/l_deltaPix 
    if ratio <0.98:
        img_show = zoom(img_org, ratio)
    else:
        img_show  = img_org
    if len(img_show)/2 == int(len(img_show)/2):
        img_show = zoom(img_show, len(img_show)/(len(img_show)+1))
    image_list[i] = img_show
    
size = np.min([len(image_list[i]) for i in range(len(image_list)) ])
for i, image in enumerate(image_list):
    # print(run_filt_list[i], image.shape, ':')
    # plt_fits(image)
    ct = int((len(image) -  size)/2)
    if ct == 0:
        continue
    else:
        image = image[ct:-ct, ct:-ct]
    # image[image<0] = 1.e-8
    # image = 2.5*np.log10(image)
    image_list[i] = image

zp_dict = {}
for i in range(len(run_filt_list)):
    mag_correct = 0
    if filt == 'F444W':
        if fit_run.fitting_specify_class.data_process_class.target_pos[0] < 5000: #module A
            correct = 0.44157708 / 0.343
            mag_correct = +2.5*np.log10(correct)
        if fit_run.fitting_specify_class.data_process_class.target_pos[0] > 5000: #module B
            correct = 0.3899884 / 0.335
            mag_correct = +2.5*np.log10(correct)
    elif filt == 'F410M':
        if fit_run.fitting_specify_class.data_process_class.target_pos[0] < 5000:
            correct = 0.9355298 / 0.832
            mag_correct = +2.5*np.log10(correct)
        if fit_run.fitting_specify_class.data_process_class.target_pos[0] > 5000:
            correct = 0.9272488 / 0.811
            mag_correct = +2.5*np.log10(correct)
    filt = run_filt_list[i]
    fit_run = fit_run_dict[filt]
    # print(filt, mag_correct)
    zp_dict[run_filt_list[i]] = fit_run.zp + mag_correct

#%%
#     # host_residual_list = fit_run.
target_id, z = load_info(idx)
smass_image = np.zeros_like(image)
sfr_image = np.zeros_like(image)
age_image = np.zeros_like(image)
AV_image = np.zeros_like(image)
import pickle
from matplotlib.colors import LogNorm      
sed_2d_info = pickle.load(open('sed_2d_info.pkl','rb'))
f = open("sed_2d_result.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
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

# norm = None  
# norm = LogNorm(vmin=4.5, vmax=8)#np.max(img[~np.isnan(img)]))
# plt.imshow(smass_image, norm=norm, origin='lower' ) 
# plt.colorbar()
# plt.show()

# norm = LogNorm(vmin=0.003, vmax=0.1)#np.max(img[~np.isnan(img)]))
# plt.imshow(sfr_image, norm=norm, origin='lower' ) 
# plt.colorbar()
# plt.show()


# norm = LogNorm(vmin=0.002, vmax=3)#np.max(img[~np.isnan(img)]))
# plt.imshow(age_image, norm=norm, origin='lower' ) 
# plt.colorbar()
# plt.show()


# # norm = LogNorm(vmin=0.001, vmax=3)#np.max(img[~np.isnan(img)]))
# norm = None
# plt.imshow(AV_image, norm=norm, origin='lower' ) 
# plt.colorbar()
# plt.show()


# import matplotlib
# cmap_r = matplotlib.cm.get_cmap('RdBu_r')

#%%
         
image_list = pickle.load(open('color_image'+'.pkl','rb'))   #use_filt = ['F444W',  'F277W', 'F150W']

from astropy.visualization import make_lupton_rgb
rgb_default = make_lupton_rgb(image_list[0][1:-1,1:-1], image_list[1][1:-1,1:-1], 
                              image_list[2][1:-1,1:-1], Q=8, stretch=0.2)
# plt.imshow(rgb_default, origin='lower')
# plt.show()

from galight.tools.measure_tools import mask_obj
from photutils import EllipticalAperture
aper = EllipticalAperture([len(smass_image)/2, len(smass_image)/2], 45, 45, theta=0)
mask1 = 1- mask_obj(smass_image, [aper])[0]
aper = EllipticalAperture([len(smass_image)/2, len(smass_image)/2], 4, 4, theta=0)
mask2 = mask_obj(smass_image, [aper])[0]
mask = mask1 * mask2

lowmass_mask = (smass_image>6.) * mask

mask[mask==0] = np.nan
lowmass_mask[lowmass_mask==0] = np.nan

import matplotlib, copy
colorname = None #'gist_heat'
my_cmap = copy.copy(matplotlib.cm.get_cmap('RdBu_r')) # copy the default cmap
my_cmap.set_bad('white')

_row = 1
fig, (axs) = plt.subplots(_row, 4, figsize=(15, 3))
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
# norm = LogNorm(vmin=4.5, vmax=8)#np.max(img[~np.isnan(img)]))


im_0 = axs[0].imshow(rgb_default , norm=None, origin='lower', vmin=6, vmax=8.1) 
axs[0].text(5,90,'F444W+F277W+F150W', c = 'white', fontsize = 16) 
# ticks = [5, 6, 7, 8]
cbar = fig.colorbar(im_0, ax=axs[0],pad=0.01, shrink=0.95, #orientation="horizontal", 
                  aspect=15)#, ticks=ticks)
cbar.set_label('M$_*$ (logM$_{\odot}$)', fontsize=20) #, rotation=270)
cbar.ax.tick_params(labelsize=20)
cbar.remove()

from galight.tools.plot_tools import scale_bar
scale_bar(axs[0], len(rgb_default), dist=1/deltaPix, text='   1"~3.75kpc', color='white')

im_1 = axs[1].imshow(smass_image * mask, norm=None, origin='lower', vmin=6, vmax=8.1, cmap = my_cmap) 
# ticks = [5, 6, 7, 8]
cbar = fig.colorbar(im_1, ax=axs[1],pad=0.01, shrink=0.95, #orientation="horizontal", 
                  aspect=15)#, ticks=ticks)
cbar.set_label('M$_*$ (logM$_{\odot}$)', fontsize=20) #, rotation=270)
cbar.ax.tick_params(labelsize=20)

# norm = LogNorm(vmin=0.001, vmax=0.01)#np.max(img[~np.isnan(img)]))
im_2 = axs[2].imshow(sfr_image * mask, norm=None, origin='lower', cmap = my_cmap, vmin = -4, vmax =-2.1)  
# ticks = [2.e-3, 6.e-3]
cbar = fig.colorbar(im_2, ax=axs[2],pad=0.01, shrink=0.95,  #orientation="horizontal", 
                  aspect=15) #, ticks=ticks)
cbar.set_label('SFR (logM$_{\odot}$/yr)', fontsize=20) #, rotation=270)
cbar.ax.tick_params(labelsize=20)

# im_4 = axs[2].imshow( (sfr_image - smass_image) * mask, norm=None, 
#                      vmin = -10.5, vmax = -8.5,
#                      origin='lower', cmap = my_cmap)  
# cbar = fig.colorbar(im_4, ax=axs[2],pad=0.01, shrink=0.95,  orientation="horizontal", 
#                   aspect=15) #, ticks=ticks)
# cbar.set_label('sSFR (logM$_{\odot}$/yr)', fontsize=20) #, rotation=270)
# cbar.ax.tick_params(labelsize=20)

# norm = LogNorm(vmin=0.5, vmax=3)#np.max(img[~np.isnan(img)]))
im_2 = axs[3].imshow(age_image * mask , norm=None, vmin=0.8, vmax=2.6, origin='lower', cmap = my_cmap) 
cbar = fig.colorbar(im_2, ax=axs[3],pad=0.01, shrink=0.95,  #orientation="horizontal", 
                   aspect=15) #, ticks=ticks)
cbar.set_label('age (Gyr)', fontsize=20) #, rotation=270)
cbar.ax.tick_params(labelsize=20)

# norm = None
# im_3 = axs[4].imshow(AV_image * mask * lowmass_mask, norm=norm, origin='lower', vmax = 4, cmap = my_cmap) 
# cbar = fig.colorbar(im_3, ax=axs[4],pad=0.01, shrink=0.95,  orientation="horizontal", 
#                   aspect=15) #, ticks=ticks)
# cbar.set_label('AV', fontsize=20) #, rotation=270)
# cbar.ax.tick_params(labelsize=20)
# ticks= np.array([1.e-4, 1.e-3, 1.e-2,1.e-1,0, 10])
# # f.colorbar(im_0, ax=axs[0], shrink=0.48, pad=0.01,  orientation="horizontal", 
# #                   aspect=15, ticks=ticks)
for _i in range(4):
    axs[_i].axis('off')
    axs[_i].axes.xaxis.set_visible(False)
    axs[_i].axes.yaxis.set_visible(False)
plt.tight_layout()
plt.subplots_adjust(wspace=-0.02, hspace=0)
plt.savefig('outcomes/SED_map.pdf')
plt.show()  

#%%
target_id, z = load_info(idx = 1)
z = float(z)

deltaPix = fit_run_dict['F356W'].fitting_specify_class.deltaPix
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
kpc_per_pixel = scale_relation * deltaPix  #kpc/pixel

sfr = 10**sfr_image  #from log(M_sun) to M_sun

def cal_sfr_dens(img, annuli = [], if_plot=False, sum_way = 'sum'):
    aper_1 = EllipticalAperture([len(img)/2, len(img)/2], 
                                annuli[0]/kpc_per_pixel, annuli[0]/kpc_per_pixel, theta=0)
    mask_1 = mask_obj(img, [aper_1])[0]  #1kpc aperture
    aper_2 = EllipticalAperture([len(img)/2, len(img)/2], 
                                annuli[1]/kpc_per_pixel, annuli[1]/kpc_per_pixel, theta=0)
    mask_2 = mask_obj(img, [aper_2])[0]  #1kpc aperture
    mask = mask_1 * (1-mask_2)
    if sum_way == 'sum':
        sum_img = np.sum(img* mask )
        img_dens = sum_img / ( np.pi *  (annuli[1] **2 - annuli[0]**2))
    elif sum_way == 'ave':
        img_dens = np.mean(img[mask!=0])
    mask[mask==0] = np.nan
    if if_plot==True:
        plt.imshow(np.log10(img)* mask, cmap = my_cmap, vmin = -4, vmax =-2.1)
        plt.show()
    return img_dens

# annuli = [1, 2.7]  #kpc
annuli = [2.7, 8]  #kpc
print( round(cal_sfr_dens(sfr, annuli=annuli, if_plot=True),3 ), 'M_sun/yr/kpc^2' )
print( round(cal_sfr_dens(age_image, annuli=annuli, if_plot=False,
                          sum_way = 'ave'),3 ), 'Gyr' )
    