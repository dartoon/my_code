#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 14:02:14 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle

import sys
sys.path.insert(0,'/Users/Dartoon/Astro/Projects/my_code/projects/2022_JWST_QSOz6/model_z6_data_id0/')

from target_info import target_info
from matplotlib.colors import LogNorm
filters = ['F356W']
# filters = ['F150W']
import copy, matplotlib

# ID_list, data_list, host_list = [], [], []

# idx = 1
data_dict = {}
host_dict = {}
redshift_dict = {}

# NOs = 10
# for idx in range(NOs):
    
# idxs = [0, 1, 2, 5, 6, 7, 9] # F356W detection
# idxs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]  #ALL
idxs = [5, 6, 7, 9]  # F150W detection
for idx in idxs:
    info = target_info[str(idx)]
    target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
    run_folder = '../model_z6_data_id{0}/stage3_all/'.format(idx)  # !!!
    # print('work on', target_id, idx)
    z_str = str(z)
    top_psf_id = 0
    count = 0
    fit_run_list = []
    filt = filters[count]
    if filt == 'F150W' :
        cmap = 'inferno'
    else:
        cmap = 'gist_heat'
    my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap)) # copy the default cmap
    my_cmap.set_bad('black')
    
    fit_files = glob.glob(
        run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))  # +\
    if idx == 1:
        fit_files = glob.glob(
            run_folder+'*fit_material*/fit_run_fixn1__idx{0}_{1}_*.pkl'.format(idx, filt))  # +\
    fit_files.sort()
    for i in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[i], 'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort()  
    fit_run = fit_run_list[sort_Chisq[top_psf_id]]
    data_dict[target_id] = fit_run.flux_2d_out['data']
    host_dict[target_id] = fit_run.flux_2d_out['data-point source']
    redshift_dict[target_id] = z

#%%
l = 0.5/0.031
def scale_bar(ax, d, dist=1/0.13, text='1"', color='black', flipped=False, fontsize=20):
        p0 = d / 15.
        ax.plot([p0, p0 + dist], [p0, p0], linewidth=2, color=color)
        ax.text(p0 + dist / 2., p0 + 0.03 * d, text, fontsize=fontsize, color=color, ha='center')

IDs = list(data_dict.keys())

if filt == 'F150W':
    my_cmap = copy.copy(matplotlib.cm.get_cmap('inferno')) # copy the default cmap
else:
    my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')
NOs = len(idxs)
fig, axes = plt.subplots(2, NOs, figsize=(NOs*2, 4))
for j in range(NOs):
    for i in range(2):
        ax = axes[i, j]
        if i == 0:
            frame_size = len(data_dict[IDs[j]])
            ax.imshow(data_dict[IDs[j]],origin='lower', cmap=my_cmap, norm=LogNorm(vmax = 20, vmin = 1.e-4))  # Change cmap and data as needed
            ax.text(frame_size*0.02, frame_size*0.9, IDs[j],fontsize=15, weight='bold', color='white')
            if j == 0:
                ax.text(frame_size*0.65, frame_size*0.1, filt,fontsize=13, weight='bold', color='white', 
                        bbox=dict(facecolor='black', alpha=0.5, linewidth=0) )
                scale_bar(ax, frame_size, dist=l, text='0.5"', color = 'white', fontsize=15)
        if i == 1:
            frame_size = len(host_dict[IDs[j]])
            ax.imshow(host_dict[IDs[j]],origin='lower', cmap=my_cmap, norm=LogNorm(vmax = 20, vmin = 1.e-4))  # Change cmap and data as needed
            ax.text(frame_size*0.05, frame_size*0.88, 'z = {0}'.format(redshift_dict[IDs[j]]),fontsize=15, weight='bold', color='white')
            if j == 0:
                if filt != 'F150W':
                    ax.text(frame_size*0.05, frame_size*0.07, 'data$-$PS = host galaxy',fontsize=12, weight='bold', color='white',
                            bbox=dict(facecolor='black', alpha=0.5, linewidth=0))
        ax.axis('off')
plt.tight_layout()
plt.subplots_adjust(wspace=-0.1, hspace=0.03)  # Adjust the values as needed
# plt.savefig('/Users/Dartoon/Downloads/host_detection_{0}.pdf'.format(filt))
plt.show()

#%%
# 0.65646

