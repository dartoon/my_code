#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:50:02 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import glob
from matplotlib.colors import LogNorm
from functions_for_result import load_info,load_prop

#%%
filts = ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']
image_list = []
fit_file_folder ='fit_material'

root_folder = '../*/*' 

# for idx in [16,20,21]:
idx = [16,20,21][2]
target_id, _ = load_info(idx)
fit_run_dict = load_prop(idx, root_folder = root_folder, prop_name='fit_run')
print(target_id)
for filt in filts:
    fit_run_list = []
    # fit_files = glob.glob(fit_file_folder+'/fit_notrunyet_{0}_idx{1}_psf0.pkl'.format(filt,idx))
    # fit_files.sort()
    f = open("./material/target_info.txt","r")
    string_1 = f.read()
    lines_ = string_1.split('\n')   # Split in to \n
    spec_z = [lines_[i].split(' ')[3] for i in range(len(lines_)) if lines_[i].split(' ')[0] == target_id][0]
    photo_z = [lines_[i].split(' ')[4] for i in range(len(lines_)) if lines_[i].split(' ')[0] == target_id][0]
    z = float(spec_z)
    if z >0:
        zinfo = 'Zspec'+str(z)
    elif z <0:
        zinfo = 'Zphot'+str(photo_z)
    fit_run = fit_run_dict[filt]
    fit_run.plot_final_qso_fit(target_ID = target_id+'_'+filt)
    print(filt,"Host mag",round(fit_run.final_result_galaxy[0]['magnitude'],3), "AGN mag",round(fit_run.final_result_ps[0]['magnitude'],3))
    image_list.append(fit_run.fitting_specify_class.data_process_class.target_stamp)


# #%%     
# import copy, matplotlib
# my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
# my_cmap.set_bad('black')
# fig, axs = plt.subplots(len(image_list), figsize=(5,25))
# for i in range(len(image_list)):
#     norm = LogNorm(vmin = 0.001, vmax =image_list[i].max() )
#     axs[i].imshow(image_list[i], norm=norm, origin='lower',cmap = my_cmap) 
#     axs[i].set_ylabel(filts[i],fontsize=20)
#     axs[i].tick_params(labelsize=15)
# check_name = target_id
# fig.suptitle('{0}'.format(check_name),fontsize=30, y=0.99)
# fig.tight_layout()
# plt.show()
