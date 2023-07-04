#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 13:23:13 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt


import sys
sys.path.insert(0, '../model_JWST_stage3_Masafusa/')
from functions_for_result import esti_smass, load_prop, load_info
from scipy.ndimage import zoom

# ID, mags, z = 'idx0', 
idx = 21
root_folder = '../model_JWST_stage3_Masafusa/'  #Include HST
fit_run_dict = load_prop(idx, root_folder = root_folder, prop_name='fit_run')
filt_list = list(fit_run_dict.keys())
# host_residual_list = []
deltaPix_list = []
# for filt in filt_list:
#     fit_run = fit_run_dict[filt]
#     # host_residual_list.append(fit_run.flux_2d_out['data-Point Source'])
#     deltaPix_list.append(fit_run.fitting_specify_class.deltaPix)
# deltaPix_list = np.array(deltaPix_list)
#%%
l_band = 'F356W'
fit_run = fit_run_dict[l_band]
l_deltaPix = fit_run.fitting_specify_class.deltaPix
# pos = fit_run.final_result_ps[0]['ra_image'], fit_run.final_result_ps[0]['dec_image'] 

from galight.tools.astro_tools import plt_fits

run_filt_list = [l_band] + [filt_list[i] for i in range(len(filt_list)) if filt_list[i] != l_band]


#%% #!!! Consider the size and the alignment issue.
size = 32  #For LW it is 56*4 and for SW it is 56*2
image_list = [None] * len(run_filt_list)
ratio_list = []
for i, filt in enumerate(run_filt_list[:]):
    fit_run = fit_run_dict[filt]
    img_org = fit_run.flux_2d_out['data'] #- np.sum(fit_run.image_host_list[1:],axis=0 )
    if len(fit_run.image_host_list) == 3:
        img_org = img_org - fit_run.image_host_list[1]
    deltaPix = fit_run.fitting_specify_class.deltaPix
    print(filt)
    #Check the host ratio is all >98%. So, no quasar subtraction:
    # print(np.sum(fit_run.image_ps_list[0]) / (np.sum(fit_run.image_ps_list[0])+ np.sum(fit_run.image_host_list[0]) ))
    # plt_fits(img_org)
    ratio = fit_run.fitting_specify_class.deltaPix/l_deltaPix 
    ratio_list.append(ratio)
    if ratio>0.8:
        pos = np.where(img_org.max() == img_org)
        pos = [pos[0][0], pos[1][0]]
        cut_size = size *2
        x_bot = pos[0]
        y_bot = pos[1]
        
    elif ratio<0.8:
        pos = np.where(img_org[30:-30,30:-30].max() == img_org)
        pos = [pos[0][0], pos[1][0]]
        cut_size = size *4
        x_bot = pos[0] + 4
        y_bot = pos[1] + 4
    ct = int(cut_size/2)
    plt_fits(img_org)
    img_org = img_org[x_bot-ct:x_bot+ct,y_bot-ct:y_bot+ct]
    # img_org = img_org[ct:-ct, ct:-ct]
    # if ratio>0.8:
    # plt_fits(img_org)
    print(pos)
    if ratio <0.98:
        img_show = zoom(img_org, 0.5)
    else:
        img_show  = img_org
        
    img_show = zoom(img_show, 0.5)
    
    img_show = img_show/np.sum(img_show) * np.sum(img_org)
    plt_fits(img_show)
    image_list[i] = img_show
    print(filt,'finish', img_show.shape)
    

#%%
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
    print(filt, mag_correct)
    zp_dict[run_filt_list[i]] = fit_run.zp + mag_correct

#%%
#     # host_residual_list = fit_run.
target_id, z = load_info(idx)
# # #%%
sed_image = np.zeros_like(image_list[0])
# write_file = open('sed_2d_info.txt','w') 
sed_2d_info = []
# write_file.write("Test of input\n")
count = 1
for i in range(len(sed_image[0])):
    for j in range(len(sed_image[1])):
        # image[i,j]
        mag_result = {}
        for k in range(len(run_filt_list)):
            filt = run_filt_list[k]
            flux = image_list[k][i,j]
            if flux>0:
                mag = -2.5*np.log10(flux) + zp_dict[filt]
                mag_result[filt] = mag
        sed_2d_info.append([i, j, mag_result])
import pickle      
pickle.dump(sed_2d_info, open('sed_2d_info_bin2.pkl', 'wb'))

pickle.dump(image_list, open('colorimage_bin2_{0}.pkl'.format(target_id), 'wb'))


#%%

images = []
filts = []
zp_list = []
for i in [-1,-3,3]:  #['F356W', 'F200W', 'F115W', 'F150W', 'F277W', 'F410M', 'F444W']
    filts.append(run_filt_list[i])
    zp_list.append(zp_dict[run_filt_list[i]])
    images.append(image_list[i])

    
from galight.tools.astro_tools import plt_fits_color
# import pickle
# pickle.dump(image_list, open('color_image_quasar'+'.pkl', 'wb'))  
images = [images[i] * 10 ** (-0.4*(zp_list[i]-zp_list[0])) for i in range(len(images)) ]
for i in range(len(images)):
    plt_fits(images[i], vmin=0.001, vmax=2.5)

plt_fits_color(images, Q=7, stretch=0.3)

#%%
from matplotlib.colors import LogNorm
import copy, matplotlib
from astropy.visualization import make_lupton_rgb
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')
def scale_bar(ax, d, dist=1/0.13, text='1"', text2=None, color='black', flipped=False, fontsize=20):
    p0 = d / 7.
    ax.plot([p0, p0 + dist], [p0, p0], linewidth=2, color=color)
    ax.text(p0 + dist / 2., p0 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')
    if text2 is not None:
        ax.text(p0 + dist / 2., p0 - 0.08 * d, text2, fontsize=fontsize, color=color, ha='center')

plot_filts = ['F150W', 'F200W', 'F277W', 'F356W']

# image_list_ct = [fit_run_dict[filt].flux_2d_out['data'] for filt in plot_filts]

fig, axs = plt.subplots(len(plot_filts)+1, figsize=(5,18))
for i in range(len(plot_filts)):
    filt = plot_filts[i]
    image_list_ct = fit_run_dict[filt].flux_2d_out['data']
    norm = LogNorm(vmin = np.std(image_list_ct[:,:2])*0.6, vmax =image_list_ct.max()/1.5 )
    axs[i].imshow(image_list_ct, norm=norm, origin='lower',cmap = my_cmap) 
    axs[i].set_ylabel(filt,fontsize=20)
    axs[i].tick_params(labelsize=15)
    pixscale = fit_run_dict[filt].fitting_specify_class.deltaPix
    axs[i].set_xticks(np.arange(0,len(image_list_ct), 1/pixscale))
    axs[i].set_yticks(np.arange(0,len(image_list_ct), 1/pixscale))
    axs[i].set_xticklabels(np.round(np.arange(0,len(image_list_ct), 1/pixscale)*pixscale,2))
    axs[i].set_yticklabels(np.round(np.arange(0,len(image_list_ct), 1/pixscale)*pixscale,2))
    # if i == 0:
    #     scale_bar(axs[i], len(image_list_ct[i]), dist=0.5/deltaPix, text='0.5"', text2 ='{0:.2f}kpc'.format(scale), color = 'white')
rgb_default = make_lupton_rgb(images[0], images[1], images[2], Q=7, stretch=0.3)

filters = [['F356W', 'F200W', 'F115W', 'F150W', 'F277W', 'F410M', 'F444W'][i] for i in [-1, -3, 3]]
use_filt = filters[0]+ '+'+ filters[1]+'+' + filters[2]

# fig, ax = plt.subplots()
axs[-1].imshow(rgb_default, origin='lower')
# plt.text(1,80,use_filt,fontsize=20, color = 'white')
pixel_s = 0.03*2
axs[-1].set_xticks(np.arange(0,len(images[0]), 1/pixel_s))
axs[-1].set_yticks(np.arange(0,len(images[0]), 1/pixel_s))
axs[-1].set_xticklabels(np.arange(0,len(images[0]), 1/pixel_s)*pixel_s)
axs[-1].set_yticklabels(np.arange(0,len(images[0]), 1/pixel_s)*pixel_s)
plt.text(1,27,use_filt,fontsize=20, color = 'white')
axs[-1].tick_params(labelsize=15)

target_id = target_id.replace('aegis_', 'AEGIS ')
fig.suptitle('{0}'.format(target_id),fontsize=35)
fig.tight_layout()
fig.savefig('/Users/Dartoon/Downloads/{0}_filt_color.pdf'.format(target_id))
plt.show()


# #%%
# # sed_2d_info = pickle.load(open('sed_2d_info.pkl','rb'))
# # count = 100
# # mag_dict = sed_2d_info[count][2]
# # esti_smass(ID = '202208'+str(count), mags_dict = mag_dict, z = z, flag = 1, if_run_gsf=True)
# import glob
# import os # Delete xfile.txt
# # folder = 'esti_smass/202208'+str(count)
# folder = 'esti_smass/2022081'
# spec_file = glob.glob(folder+'/gsf_spec_*.fits')[0]
# hdul_spec = pyfits.open(spec_file)
# info_spec = hdul_spec[1].header

# write_file = open(folder + '/gsf_spec_header.txt','w') 
# write_file.write(str(info_spec))
# write_file.close()
# rm_file = glob.glob(folder+'/gsf_spec_*.fits') + glob.glob(folder + '/SPEC*png') + glob.glob(folder+'/*asdf') 
# for file in rm_file:
#     os.remove(file)
# steller_file = glob.glob(folder+'/SFH_*.fits')[0]
# hdul = pyfits.open(steller_file)
# info = hdul[0].header 
# smass = info['Mstel_50']
# sfr = info['SFR_50']
# m_age = info['T_MW_50']
# l_age = info['T_LW_50']
# AV = info['AV_50']
# filename = 'sed_2d_result.txt'
# if_file = glob.glob(filename)
# if if_file == []:
#     write_file =  open(filename,'w')
#     write_file.write("count_i, smass, sfr, m_age, l_age, AV \n")
# else:
#     write_file =  open(filename,'r+') 
#     write_file.read()
# write_file.write("{0} {1} {2} {3} {4} {5}".format(count, smass, sfr, m_age, l_age, AV))
# write_file.write("\n")
# write_file.close()

# spec_file = glob.glob('esti_smass/2022081/gsf_spec_*.fits')[0]
# hdul_spec = pyfits.open(spec_file)
# info_spec = hdul_spec[1].header

# JWST_smass.append( float(info['Mstel_50']) )