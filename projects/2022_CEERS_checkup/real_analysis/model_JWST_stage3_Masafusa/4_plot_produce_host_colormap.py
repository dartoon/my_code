#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 20:30:48 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import os
import glob

remove_id = [24, 55]
# use_filt = ['F444W', 'F410M', 'F356W', 'F277W', 'F200W', 'F150W', 'F115W', 'F160W' ,'F125W' ,'F814W', 'F606W']

# use_filt = ['F356W',  'F277W', 'F150W']
use_filt = ['F444W',  'F277W', 'F150W']   
fit_run_list = []

for idx in [51]:  #z_spec > 1.6
    if idx in remove_id:
        continue
    else:
        files = glob.glob('../*/fit_material/data_process_idx{0}_*_psf*.pkl'.format(idx))
        files.sort()
        
        _collect_info = []
        for i in range(len(files)):
            _file = files[i]
            idx_info = _file.split('idx')[1].split('_')[0]
            filt_info = _file.split('_psf')[0].split('_')[-1]
            this_info = [idx_info, filt_info]
            if this_info not in _collect_info:
                _collect_info.append(this_info)
        
        filters = [_collect_info[i][1] for i in range(len(_collect_info))]
        if 'F814W' in filters:
            filters = ['F814W'] + [filters[i] for i in range(len(filters)) if filters[i]!='F814W' ]
        if 'F606W' in filters:
            filters = ['F606W'] + [filters[i] for i in range(len(filters)) if filters[i]!='F606W' ]
        
        f = open("../model_JWST_stage3_Masafusa/target_idx_info.txt","r")
        string = f.read()
        lines = string.split('\n')   # Split in to \n
        target_id = [lines[i].split(' ')[1] for i in range(len(lines)) if lines[i].split(' ')[0] == str(idx)][0]
        
        for filt in use_filt:
            _fit_run_list = []
            idx = idx_info
            # idx, filt= item
            fit_files = glob.glob('../model_JWST_stage3_Masafusa/fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))+\
                        glob.glob('../model_HST_material/fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))
            fit_files.sort()
            warn_strs = ['F115W_psf6', 'F150W_psf7', 'F277W_psf2']
            for warn_str in warn_strs:
                fit_files = [fit_files[i] for i in range(len(fit_files)) if warn_str not in fit_files[i]]
            
            for i in range(len(fit_files)):
                _fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
            chisqs = np.array([_fit_run_list[i].reduced_Chisq for i in range(len(_fit_run_list))])
            sort_Chisq = chisqs.argsort()  
            print('idx', idx, filt, "Total PSF NO.", len(sort_Chisq))
            
            if 'HST_material' in fit_files[0]:
                weight = np.zeros(len(chisqs))
                weight[sort_Chisq[0]] = 1
            elif 'JWST_stage3' in fit_files[0]:
                count_n = 5
                Chisq_best = chisqs[sort_Chisq[0]]
                Chisq_last= chisqs[sort_Chisq[count_n-1]]
                inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
                weight = np.zeros(len(chisqs))
                for i in sort_Chisq[:count_n]:
                    weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
            fit_run = _fit_run_list[sort_Chisq[0]]
            fit_run_list.append(fit_run)
            # fit_run.plot_final_qso_fit(target_ID = target_id+'-'+filt)
            prop_name = 'n_sersic'
            all_values = [_fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(_fit_run_list))]
            weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
            rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
            host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
            AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
            ratio = host_flux/(host_flux+AGN_flux)
            print(prop_name, round(weighted_value,2), '+-', round(rms_value,2))
            print('Chisqs top 2', round(chisqs[sort_Chisq[0]],2), round(chisqs[sort_Chisq[1]],2))
            print_s =filt +' ratio: ' + str(round(ratio,2)) + "\n\n\n"
            print(fit_files[sort_Chisq[0]])
            print(print_s)
            # hold = input(print_s)

zp_list = [fit_run_list[i].zp for i in range(3)]

#%%
from scipy.ndimage import zoom
shift_center = True
l_idx = 0
deltaPix_list = np.array([fit_run_list[i].fitting_specify_class.deltaPix for i in range(len(fit_run_list))])
ratio_list = deltaPix_list/np.max(deltaPix_list)
image_list = []
from galight.tools.astro_tools import plt_fits
pos = fit_run.final_result_ps[0]['ra_image'], fit_run.final_result_ps[0]['dec_image'] 
run_idx_list = [l_idx] + [i for i in range(len(use_filt)) if i != l_idx]

image_list = [None] * len(use_filt)
for i in run_idx_list:
    fit_run = fit_run_list[i]
    # img_org = fit_run.flux_2d_out['data'] 
    img_org = fit_run.flux_2d_out['data-Point Source'] 
    if len(fit_run.image_host_list) == 3:
        img_org = img_org - fit_run.image_host_list[1]
    # img_org = fit_run.flux_2d_out['model']
    if shift_center == True:
        if hasattr(fit_run, 'final_result_ps'):
            pos =  - fit_run.final_result_ps[0]['ra_image'][0]/deltaPix_list[i] ,\
                  fit_run.final_result_ps[0]['dec_image'][0]/deltaPix_list[i]
        else:
            pos =  - fit_run.final_result_galaxy[0]['center_x']/deltaPix_list[i] ,\
                  fit_run.final_result_galaxy[0]['center_y']/deltaPix_list[i]
        pos = np.int0(np.array(pos))
        if i == l_idx:
            l_pos = pos
        if i != l_idx:
            shift = pos - l_pos
            new_pos = pos + shift + int(len(img_org)/2) #!!!
            ct = int(len(img_org)/2) - np.max(shift) - 1 
            img_org =  img_org[ new_pos[1] - ct:new_pos[1] + ct+1, new_pos[0] - ct:new_pos[0] + ct +1 ]
    print(ratio_list[i])
    img_show = zoom(img_org, ratio_list[i])
    print(img_show.shape)
    if len(img_show)/2 != int(len(img_show)/2):
        img_show = zoom(img_show, len(img_show)/(len(img_show)-1))
    img_show = img_show/np.sum(img_show)*np.sum(img_org)
    print(use_filt[i])
    plt_fits(img_show)
    print(i)
    image_list[i] = img_show
#%%
size = np.min([len(image_list[i]) for i in range(len(image_list)) ])
for i, image in enumerate(image_list):
    ct = int((len(image) -  size)/2)
    if ct == 0:
        continue
    else:
        image = image[ct:-ct, ct:-ct]
    # image[image<0] = 1.e-8
    # image = 2.5*np.log10(image)
    image_list[i] = image
    
image_list = [image_list[i] * 10 ** (-0.4*(zp_list[i]-zp_list[0])) for i in range(3) ]
    
from galight.tools.astro_tools import plt_fits_color
# pickle.dump(image_list, open('color_image_quasar'+'.pkl', 'wb'))  
plt_fits_color(image_list, Q=7, stretch=0.3)

# from astropy.visualization import make_lupton_rgb
# rgb_default = make_lupton_rgb(image_list[0], image_list[1], image_list[2], Q=8, stretch=0.2)
# plt.imshow(rgb_default, origin='lower')
# plt.show()
