#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 21:31:32 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from functions_for_result import esti_smass, load_prop, load_info
from scipy.ndimage import zoom
import copy, matplotlib
from matplotlib.colors import LogNorm
import glob
import pickle 
from matplotlib.ticker import AutoMinorLocator
from galight.tools.measure_tools import SB_profile

idx = 35

fit_files = glob.glob('./*fit_material/fit_run_idx{0}_{1}_CombPsfsNO_8_*.pkl'.format(idx, '*'))
filt_list = ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
fig, (axs) = plt.subplots(2, 4, figsize=(15, 8) )
# for i in range(len(imgs)):
    
from galight.tools.measure_tools import measure_FWHM
# for filt
for i_filt, filt in enumerate(filt_list):
    fit_file = [fit_files[i] for i in range(len(filt_list)) if filt in fit_files[i]][0]
    fit_run = pickle.load(open(fit_file,'rb'))
    PSF_list = fit_run.fitting_specify_class.data_process_class.PSF_list
    prf_name_list = ['PSF{0}'.format(i) for i in range(len(PSF_list))]
    webbPSF = pyfits.getdata('webbPSFs/PSF_{0}.fits'.format(filt))
    PSF_list.append(webbPSF)
    prf_name_list.append('webbpsf')
    # profiles_compare(PSF_list, prf_name_list = prf_name_list, if_annuli= False)
    _i = int(i_filt / 3)
    _j = int(i_filt % 3)
    if _i == 2 and _j ==0:
        _i = 1
        _j = 3

    prf_list = PSF_list
    prf_NO = len(prf_list)
    scale_list = [1] * len(prf_list)
    grids = 20  
    norm_pix = 3 
    if_annuli=True
    y_log=False
    x_gridspace = None
    hide_axes=False
    if len(prf_name_list)!=len(prf_list):
        raise ValueError("The profile name is not in right length")
    FWHM = np.mean( measure_FWHM(prf_list[0]))
    FWHM_webbpsf = np.mean( measure_FWHM(prf_list[-1]))
    for i in range(prf_NO):
        b_c = int(len(prf_list[i])/2)
        b_r = int(len(prf_list[i])/6)
        center = np.reshape(np.asarray(np.where(prf_list[i]== prf_list[i][b_c-b_r:b_c+b_r,b_c-b_r:b_c+b_r].max())),(2))[::-1]
        scale = scale_list[i]
        if _i == 0:
            radius = 12
        if _i > 0:
            radius = 6
        
        r_SB, r_grids = SB_profile(prf_list[i], center, radius=radius*scale,
                                    grids=grids, x_gridspace=x_gridspace,if_annuli=if_annuli)
        if isinstance(norm_pix,int) or isinstance(norm_pix,float):
            count = r_grids <= norm_pix * scale
            idx = count.sum() -1
#            print("idx:",idx)
            r_SB /= r_SB[idx]      #normalize the curves
        r_grids /= scale
        if i < prf_NO-1:
            if y_log == False:
                axs[_i][_j].plot(r_grids, r_SB, 'x-', label=prf_name_list[i])
            elif y_log == True:
                axs[_i][_j].plot(r_grids, np.log10(r_SB), 'x-', label=prf_name_list[i])
        elif i == prf_NO-1:
            if y_log == False:
                axs[_i][_j].plot(r_grids, r_SB, '-',c='skyblue' , label=prf_name_list[i], linewidth = 6, alpha = 0.9)
            elif y_log == True:
                axs[_i][_j].plot(r_grids, np.log10(r_SB), '--',c='skyblue', label=prf_name_list[i], linewidth = 6, alpha = 0.9)
            axs[_i][_j].plot(np.linspace(0,10) * 0 + FWHM, np.linspace(0,10) , '--',c='red', label='FWHM\nstar', linewidth = 3, alpha = 1)
            axs[_i][_j].plot(np.linspace(0,10) * 0 + FWHM_webbpsf, np.linspace(0,10) , '--',c='skyblue', label='FWHM\nwebbPSF', linewidth = 3, alpha = 1)
            
    axs[_i][_j].tick_params(which='both', width=2, direction='in')
    axs[_i][_j].tick_params(which='major', length=7, direction='in')
    axs[_i][_j].tick_params(which='minor', length=4, color='r', direction='in')
    axs[_i][_j].text(4.9, 2.0,filt, fontsize=25)
    if _i == 0:
        axs[_i][_j].set_xlim([1.4, 6.9])
        axs[_i][_j].set_ylim([0.04, 2.5])
    if _i == 1:
        axs[_i][_j].set_xlim([1.4, 6.9])
        axs[_i][_j].set_ylim([0.04, 2.5])
    
    if _i == 1:
        axs[_i][_j].set_xlabel("Pixels", fontsize=20)
    if _j == 0:
        axs[_i][_j].set_ylabel("Scaled SB (annuli)", fontsize=17)
    if _j > 0:
        # axs[_i][_j].axes.yaxis.set_visible(False)
        # plt.grid(which="minor")
        axs[_i][_j].set_yticklabels([])
    if x_gridspace == 'log':
        axs[_i][_j].set_xscale('log')
        # plt.xlim(1.3, ) 
    # axs[_i][_j].xaxis.set_minor_locator(minorLocator)
    # plt.grid(which="minor")
    # plt.legend(prop={'size':15},ncol=2)
    plt.tick_params(labelsize=20)
    axs[_i][_j].tick_params(labelsize=20)      
    
    # if hide_axes == True:
    #     axs[_i][_j].axes.xaxis.set_visible(False)
    #     axs[_i][_j].axes.yaxis.set_visible(False)
        
axs[0][2].legend(prop={'size':15},ncol=2, bbox_to_anchor=(0.98, 0.15))

# plt.tick_params(labelsize=20)        
axs[0][3].axes.xaxis.set_visible(False)
axs[0][3].axes.yaxis.set_visible(False)
axs[0][3].axis('off')
plt.subplots_adjust(wspace=0.0, hspace=0.2)
plt.savefig('outcomes/psfs_profile.pdf')
plt.show()
# # for i in range( 5 - len(imgs)%5 ):
# #     axs[-1][-(i+1)].axis('off')
# if savename is not None:
#     plt.savefig(savename,bbox_inches='tight')
# if if_plot == True:
#     plt.show()    
# else:
#     plt.close()
    
    #%%
