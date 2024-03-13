#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 11:31:08 2024

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")
import glob, pickle
import sys
from matplotlib.colors import LogNorm
from photutils.aperture import EllipticalAperture
import copy, matplotlib

sys.path.insert(0, '../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info
run_folder = '../material/fit_result/'

# Set the dimensions and the number of subplots
fig, axs = plt.subplots(4, 8, figsize=(16, 8))

# Adjust the space between the subplots to zero
plt.subplots_adjust(wspace=0, hspace=0.03)
# Populate the subplots with content and remove ticks
len_j = 8
for jj in range(len_j):
    idx = jj + 2
    info = target_info[str(idx)]
    target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
    fit_run_list_sp1 = []
    add_cond = '_fixn1'
    filt = 'F356W'
    count_n = 5
    all_library = True
    fit_run_list_sp1 = []
    psf_sp = 1
    if idx == 1:  #!!!
        all_library = False
    if all_library == True:
        load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx*_psfsf{2}{3}.pkl'.format(filt, idx, psf_sp, add_cond))
    else:
        psf_idx = idx
        load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx{2}*_psfsf{3}{4}.pkl'.format(filt, idx, psf_idx, psf_sp,add_cond))
    load_files.sort()
    chisqs_idx = []
    for file in load_files:
        fit_run_list_sp1.append(pickle.load(open(file,'rb')))
    chisqs = np.array([fit_run_list_sp1[i].reduced_Chisq for i in range(len(fit_run_list_sp1))])
    sort_Chisq_sp1 = chisqs.argsort()  
    weight_sp1 = np.zeros(len(chisqs))
    for i in sort_Chisq_sp1[:count_n]:
        weight_sp1[i] = 1
    psf_sp = 2
    fit_run_list_sp2 = []
    if all_library == True:
        load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx*_psfsf{2}{3}.pkl'.format(filt, idx, psf_sp, add_cond))
    else:
        psf_idx = idx
        load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx{2}*_psfsf{3}{4}.pkl'.format(filt, idx, psf_idx, psf_sp,add_cond))
    load_files.sort()
    chisqs_idx = []
    for file in load_files:
        fit_run_list_sp2.append(pickle.load(open(file,'rb')))
    chisqs = np.array([fit_run_list_sp2[i].reduced_Chisq for i in range(len(fit_run_list_sp2))])
    sort_Chisq_sp2 = chisqs.argsort()  
    weight_sp2 = np.zeros(len(chisqs))
    for i in sort_Chisq_sp2[:count_n]:
        weight_sp2[i] = 1
    weight = np.concatenate([weight_sp1, weight_sp2])
    fit_run_list = fit_run_list_sp1 + fit_run_list_sp2    
    
    fit_run = fit_run_list_sp1[sort_Chisq_sp1[0]]
    if fit_run.reduced_Chisq > fit_run_list_sp2[sort_Chisq_sp2[0]].reduced_Chisq:
        fit_run = fit_run_list_sp2[sort_Chisq_sp2[0]]  #Save Top result
    
    vmin = 0.0001
    vmax = 3
    norm = LogNorm(vmin=vmin, vmax=vmax)#np.max(img[~np.isnan(img)]))
    my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
    my_cmap.set_bad('black')
    
    
    axs[0, jj].set_title(target_id, fontsize=17, fontweight="bold")  # Add title
    im0 = axs[0, jj].imshow(fit_run.flux_2d_out['data'], 
                      norm=norm, origin='lower',cmap = my_cmap) 
    im0 = axs[1, jj].imshow(fit_run.flux_2d_out['data-point source'], 
                      norm=norm, origin='lower',cmap = my_cmap) 
    im1 = axs[2, jj].imshow(fit_run.image_host_list[0],
                      norm=norm, origin='lower',cmap = my_cmap) 
    
    used_PSF_list = [fit_run_list[i].fitting_specify_class.data_process_class.PSF_list[0] for i in range(len(fit_run_list)) if weight[i] != 0]
    PSF_noise_map = np.std(used_PSF_list, axis=0) * fit_run.final_result_ps[0]['flux_within_frame']
    total_noise = np.sqrt(PSF_noise_map[20:-20, 20:-20] ** 2 +  fit_run.fitting_specify_class.data_process_class.noise_map**2)
    fit_run.cal_astrometry()
    x = fit_run.final_result_galaxy[0]['position_xy'][0] + len(fit_run.image_host_list[0])/2
    y = fit_run.final_result_galaxy[0]['position_xy'][1] + len(fit_run.image_host_list[0])/2
    a = fit_run.final_result_galaxy[0]['R_sersic']/fit_run.fitting_specify_class.deltaPix * 4
    if a > len(fit_run.image_host_list[0])/2 * 0.8:
        a = len(fit_run.image_host_list[0])/2 * 0.8
    b = a*fit_run.final_result_galaxy[0]['q']
    theta = - fit_run.final_result_galaxy[0]['phi_G']
    aperture = EllipticalAperture((x, y), a, b, theta=theta)
    im2 = axs[3, jj].imshow( (fit_run.flux_2d_out['data-point source'] - np.sum(fit_run.image_host_list[1:],axis = 0) ) /total_noise,
                      vmin = 0.001, vmax = 5, origin='lower')
    aperture.plot(color= 'white',
                  lw=1, label = 'comp {0}'.format(i), axes = axs[3, jj])
    cbar_ax = fig.add_axes([0.905, 0.33, 0.014, 0.55])  # Adjust these values as needed
    cbar_ax.tick_params(labelsize=14)
    fig.colorbar(im0, cax=cbar_ax)
    cbar_ax_bottom = fig.add_axes([0.905, 0.125, 0.014, 0.19])  # Adjust these values as needed
    cbar_ax_bottom.tick_params(labelsize=10)
    # # Add another colorbar to the bottom of the figure
    fig.colorbar(im2, cax=cbar_ax_bottom) #, orientation='horizontal')
    axs[0, 0].set_ylabel('Data', fontsize=15, fontweight="bold")  # Add title
    axs[1, 0].set_ylabel('Data-QSO', fontsize=15, fontweight="bold")  # Add title
    axs[2, 0].set_ylabel('Host model', fontsize=15, fontweight="bold")  # Add title
    axs[3, 0].set_ylabel('SNR', fontsize=15, fontweight="bold")  # Add title
    for ii in range(0,4):
        axs[ii, jj].set_xticks([])  # Remove x ticks
        axs[ii, jj].set_yticks([])  # Remove y ticks
# plt.savefig('/Users/Dartoon/Downloads/figure2.pdf', bbox_inches = "tight")
plt.show()

# # Show the figure
#     import copy, matplotlib
#     if cmap == 'gist_heat':
#         my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
#         my_cmap.set_bad('black')
#     else:
#         my_cmap = None
#     fig, ax = plt.subplots(figsize=figsize)
#     if norm == 'log':
#         norm = LogNorm(vmin=vmin, vmax=vmax)#np.max(img[~np.isnan(img)]))
#     else:
#         norm = norm
#     plt.imshow(img, norm=norm, origin='lower',cmap = my_cmap) 
#     if colorbar == True:
#         plt.colorbar()
#     if savename is not None:
#         plt.savefig(savename,bbox_inches='tight')
#     if hold == False:
#         plt.show()     


