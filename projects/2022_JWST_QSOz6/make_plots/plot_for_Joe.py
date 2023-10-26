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
idxs = [1]  # F150W detection
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
def scale_bar(ax, d, dist=1/0.13, text='1"', color='black', flipped=False, fontsize=20):
    if flipped:
        p0 = d - d / 15. - dist
        p1 = d / 15.
        ax.plot([p0, p0 + dist], [p1, p1], linewidth=2, color=color)
        ax.text(p0 + dist / 2., p1 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')
    else:
        p0 = d / 15.
        ax.plot([p0, p0 + dist], [p0, p0], linewidth=2, color=color)
        ax.text(p0 + dist / 2., p0 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')

my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')

from galight.tools.astro_tools import plt_fits
from matplotlib.colors import LogNorm
fig, ax = plt.subplots(figsize=None)
norm = LogNorm(vmax = 8, vmin = 1.e-4)
plt.imshow(fit_run.flux_2d_out['data-point source'], origin='lower',cmap = my_cmap, norm = norm)
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)
deltaPix = fit_run.fitting_specify_class.deltaPix

scale_bar(ax, len(fit_run.flux_2d_out['data']), dist=0.5/deltaPix, text='0.5"', color = 'w')
cb_i = plt.colorbar(#shrink=0.48, pad=0.01,  #orientation="horizontal", 
                  aspect=15)
cb_i.ax.tick_params(labelsize=15) 
plt.savefig('{0}_host.pdf'.format(target_id))
# scale_bar(ax, frame_size, dist=0.5/deltaPix, text='0.5"', color = color)
