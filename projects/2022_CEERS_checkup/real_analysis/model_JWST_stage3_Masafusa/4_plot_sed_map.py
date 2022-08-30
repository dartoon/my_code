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
root_folder = './*'  #Include HST
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
    img_org = fit_run.flux_2d_out['data-Point Source'] - np.sum(fit_run.image_host_list[1:],axis=0 )
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
    ct = int((len(image) -  size)/2)
    if ct == 0:
        continue
    else:
        image = image[ct:-ct, ct:-ct]
    # image[image<0] = 1.e-8
    # image = 2.5*np.log10(image)
    image_list[i] = image
    print(run_filt_list[i], image.shape, ':')
    plt_fits(image)

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
sed_image = np.zeros_like(image)
import pickle
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
        sed_image[_i, _j] = float(smass)
plt_fits(sed_image)
    