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
import sys
sys.path.insert(0, '../model_JWST_stage3_Masafusa/')
# ID, mags, z = 'idx0', 
# 1,2,0,51,35
idx = 16
# root_folder = '../*/*'  #Include HST
root_folder = '../model_JWST_stage3_Masafusa/'  #Not include HST
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
import pickle
from matplotlib.colors import LogNorm      
bin_info = '_bin2'
sed_2d_info = pickle.load(open('sed_2d_info{0}.pkl'.format(bin_info),'rb'))
f = open("sed_2d_result{0}.txt".format(bin_info),"r")
string = f.read()
lines = string.split('\n')   # Split in to \n
size = int(np.sqrt(len(sed_2d_info)))
smass_image = np.zeros([size,size])
sfr_image =  np.zeros([size,size])
age_image =  np.zeros([size,size])
AV_image = np.zeros([size,size])
EBV_image = np.zeros([size,size])
for ct, line in enumerate(lines[1:-1]):
    if len(line.split(' ')) < 4:
        continue
    else:
        count, smass, sfr, m_age, l_age, AV, EBV = line.split(' ')
        count = int(count)
        _i, _j = sed_2d_info[count][0], sed_2d_info[count][1]
        smass_image[_i, _j] = float(smass)    #smass in logMsun
        sfr_image[_i, _j] = float(sfr)          #logMsun/yr 
        age_image[_i, _j] = 10**float(m_age)    #logGyr to Gry
        AV_image[_i, _j] = AV    #logGyr
        EBV_image[_i, _j] = EBV    #logGyr

        
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

# norm = LogNorm(vmin=0.001, vmax=3)#np.max(img[~np.isnan(img)]))
norm = None
plt.imshow(EBV_image, norm=norm, origin='lower' ) 
plt.colorbar()
plt.show()


import matplotlib
cmap_r = matplotlib.cm.get_cmap('RdBu_r')

# #%%
# from astropy.visualization import make_lupton_rgb
# image_list_qso = pickle.load(open('color_image_quasar'+'.pkl','rb'))   #use_filt = ['F444W',  'F277W', 'F150W']
# rgb_default_qso = make_lupton_rgb(image_list_qso[0][1:-1,1:-1], image_list_qso[1][1:-1,1:-1], 
#                               image_list_qso[2][1:-1,1:-1], Q=10, stretch=0.4)
# image_list = pickle.load(open('color_image'+'.pkl','rb'))   #use_filt = ['F444W',  'F277W', 'F150W']
# rgb_default = make_lupton_rgb(image_list[0][1:-1,1:-1], image_list[1][1:-1,1:-1], 
#                               image_list[2][1:-1,1:-1], Q=8, stretch=0.2)
# plt.imshow(rgb_default, origin='lower')
# plt.show()

# #%%
# if bin_info != '':
#     dvd_info = int(bin_info[-1])
# else:
#     dvd_info = 1
    
# rad1, rad2 = 4, 8
# rad0 = 2
# target_id, z = load_info(idx = 1)
# z = float(z)
# deltaPix = fit_run_dict['F356W'].fitting_specify_class.deltaPix
# _deltaPix = deltaPix * dvd_info  #The true Pixel size after binning 
# from astropy.cosmology import FlatLambdaCDM
# cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
# scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
# kpc_per_pixel = scale_relation * _deltaPix  #kpc/pixel

# from galight.tools.measure_tools import mask_obj
# from photutils import EllipticalAperture
# aper = EllipticalAperture([len(smass_image)/2, len(smass_image)/2], 45/dvd_info, 
#                           45/dvd_info, theta=0)
# #45/dvd_info / kpc_per_pixel, outer region.
# mask1 = 1- mask_obj(smass_image, [aper])[0]
# aper = EllipticalAperture([len(smass_image)/2, len(smass_image)/2], 
#                           round(rad0/kpc_per_pixel), round(rad0/kpc_per_pixel), theta=0)
# mask2 = mask_obj(smass_image, [aper])[0]


# aper27 = EllipticalAperture([len(smass_image)/2, len(smass_image)/2], round(rad1/kpc_per_pixel), 
#                           round(rad1/kpc_per_pixel), theta=0)
# aper80 = EllipticalAperture([len(smass_image)/2, len(smass_image)/2], round(rad2/kpc_per_pixel), 
#                           round(rad2/kpc_per_pixel), theta=0)

# aper_ = EllipticalAperture([len(image_list[0])/2-1, len(image_list[0])/2-1], 
#                           round(rad0/kpc_per_pixel * 2), round(rad0/kpc_per_pixel * 2), theta=0)
# aper27_ = EllipticalAperture([len(image_list[0])/2-1, len(image_list[0])/2-1], 
#                           round(rad1/kpc_per_pixel * 2), round(rad1/kpc_per_pixel * 2), theta=0)
# aper80_ = EllipticalAperture([len(image_list[0])/2-1, len(image_list[0])/2-1], 
#                           round(rad2/kpc_per_pixel * 2), round(rad2/kpc_per_pixel * 2), theta=0)

# mask = mask1 * mask2

# lowmass_mask = (smass_image>6.) * mask

# mask[mask==0] = np.nan
# lowmass_mask[lowmass_mask==0] = np.nan

# import matplotlib, copy
# colorname = None #'gist_heat'
# my_cmap_mass = copy.copy(matplotlib.cm.get_cmap('afmhot')) # copy the default cmap
# my_cmap_mass.set_bad('white')
# my_cmap = copy.copy(matplotlib.cm.get_cmap('RdBu_r')) # copy the default cmap
# my_cmap.set_bad('white')
# my_cmap_ = copy.copy(matplotlib.cm.get_cmap('RdBu')) # copy the default cmap
# my_cmap_.set_bad('white')


# import imageio
# filename_list = []
# plt.figure(figsize=(11, 11))
# plt.title("F444W+F277W+F150W", fontsize=35)
# plt.imshow(rgb_default_qso , norm=None, origin='lower') 
# plt.axis('off')
# plt.savefig('gifs/SDSS1420A_color.png')#, dpi=1200)
# plt.close()
# plt.figure(figsize=(11, 11))
# plt.title("F444W+F277W+F150W", fontsize=35)
# plt.imshow(rgb_default , norm=None, origin='lower') 
# plt.axis('off')
# plt.savefig('gifs/SDSS1420A_color_host.png')#, dpi=1200)
# plt.close()

# filename_list = ['gifs/SDSS1420A_color.png', 'gifs/SDSS1420A_color_host.png']
# with imageio.get_writer('gifs/qso_to_host.gif', mode='I', fps=1) as writer:
#     for filename in filename_list:
#         image = imageio.imread(filename)
#         writer.append_data(image)

# #%%
# _row = 1
# fig, (axs) = plt.subplots(_row, 5, figsize=(15, 4))
# import matplotlib as mat
# mat.rcParams['font.family'] = 'STIXGeneral'
# # norm = LogNorm(vmin=4.5, vmax=8)#np.max(img[~np.isnan(img)]))

# im_0 = axs[0].imshow(rgb_default_qso , norm=None, origin='lower') 
# axs[0].text(5, 5, 'SDSS1420A',fontsize=14, 
#              weight='bold', color='white', bbox={'facecolor': 'gold', 'alpha': 0.6, 'pad': 3})

# axs[0].text(1,83,'F444W+F277W+F150W\nquasar', c = 'white', fontsize = 15) 
# # ticks = [5, 6, 7, 8]
# # cbar = fig.colorbar(im_0, ax=axs[0],pad=0.01, shrink=0.95, #orientation="horizontal", 
# #                   aspect=15)#, ticks=ticks)
# c_list = 'white'
# im_0 = axs[1].imshow(rgb_default , norm=None, origin='lower') 
# axs[1].text(1,83,'host\nF444W+F277W+F150W\nhost galaxy', c = 'white', fontsize = 15) 
# ticks = [5, 6, 7, 8]
# aper_.plot(color= c_list,
#               lw=1.0, label = 'comp {0}'.format(i), axes = axs[0], alpha =0.6)
# aper27_.plot(color= c_list,
#               lw=1.0, label = 'comp {0}'.format(i), axes = axs[0], alpha =0.6)
# aper80_.plot(color= c_list,
#               lw=1.0, label = 'comp {0}'.format(i), axes = axs[0], alpha =0.6)
# # aper_.plot(color= c_list,
# #               lw=1.0, label = 'comp {0}'.format(i), axes = axs[1], alpha =0.3)
# # aper27_.plot(color= c_list,
# #               lw=1.0, label = 'comp {0}'.format(i), axes = axs[1], alpha =0.3)
# # aper80_.plot(color= c_list,
# #               lw=1.0, label = 'comp {0}'.format(i), axes = axs[1], alpha =0.3)

# cbar = fig.colorbar(im_0, ax=axs[1],pad=0.01, shrink=0.95, #orientation="horizontal", 
#                   aspect=15)#, ticks=ticks)
# cbar.set_label('M$_*$ (logM$_{\odot}$)', fontsize=20) #, rotation=270)
# cbar.ax.tick_params(labelsize=20)
# cbar.remove()


# def scale_bar(ax, d, dist=1/0.1, text='1"', color='black', fontsize=15):
#     p0 = d - d /4. - dist/3*2
#     p1 = 3 * d / 20.
#     ax.plot([p0, p0 + dist], [p1, p1], linewidth=2, color=color)
#     ax.text(p0 + dist, p1/2 - 0.02 * d, text, fontsize=fontsize, color=color, ha='center')
        
# scale_bar(axs[0], len(rgb_default),  dist=0.24/deltaPix, text='2kpc~0.24"', color='white')

# scale_bar(axs[1], len(rgb_default),  dist=0.24/deltaPix, text='2kpc~0.24"', color='white')

# # c_list = np.random.uniform(0, 1), np.random.uniform(0, 1), np.random.uniform(0, 1)


# im_1 = axs[2].imshow(smass_image * mask, norm=None, origin='lower', vmin=7, vmax=9.8, cmap = my_cmap_mass) 
# aper27.plot(color= c_list,
#               lw=1.5, label = 'comp {0}'.format(i), axes = axs[2],alpha =0.7)
# aper80.plot(color= c_list,
#               lw=1.5, label = 'comp {0}'.format(i), axes = axs[2],alpha =0.7)
# cbar = fig.colorbar(im_1, ax=axs[2],pad=0.01, shrink=0.95, orientation="horizontal", 
#                   aspect=15)#, ticks=ticks)
# cbar.set_label('M$_*$ (logM$_{\odot}$)', fontsize=20) #, rotation=270)
# cbar.ax.tick_params(labelsize=20)

# # norm = LogNorm(vmin=0.001, vmax=0.01)#np.max(img[~np.isnan(img)]))
# im_2 = axs[3].imshow(sfr_image * mask, norm=None, origin='lower', cmap = my_cmap_, vmin = -3, vmax =-1)  
# aper27.plot(color= c_list,
#               lw=2.5, label = 'comp {0}'.format(i), axes = axs[3])
# aper80.plot(color= c_list,
#               lw=2.5, label = 'comp {0}'.format(i), axes = axs[3])

# cbar = fig.colorbar(im_2, ax=axs[3],pad=0.01, shrink=0.95,  orientation="horizontal", 
#                   aspect=15) #, ticks=ticks)
# cbar.set_label('SFR (logM$_{\odot}$/yr)', fontsize=20) #, rotation=270)
# cbar.ax.tick_params(labelsize=20)

# im_2 = axs[4].imshow(age_image * mask , norm=None, vmin=0.8, vmax=2.6, origin='lower', cmap = my_cmap) 
# aper27.plot(color= c_list,
#               lw=2.5, label = 'comp {0}'.format(i), axes = axs[4])
# aper80.plot(color= c_list,
#               lw=2.5, label = 'comp {0}'.format(i), axes = axs[4])
# cbar = fig.colorbar(im_2, ax=axs[4],pad=0.01, shrink=0.95,  orientation="horizontal", 
#                    aspect=15) #, ticks=ticks)
# cbar.set_label('age (Gyr)', fontsize=20) #, rotation=270)
# cbar.ax.tick_params(labelsize=20)

# for _i in range(5):
#     axs[_i].axis('off')
#     axs[_i].axes.xaxis.set_visible(False)
#     axs[_i].axes.yaxis.set_visible(False)
# # plt.tight_layout()
# plt.subplots_adjust(wspace=0.05, hspace=0)
# plt.savefig('outcomes/SED_map.pdf')
# plt.show()  
    