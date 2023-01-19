#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:58:31 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle

idx = 3

import sys
sys.path.insert(0, '../model_z6_data_id0/')
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# files = glob.glob('./*fit_material*sm*/data_process_idx{0}_*_psf*.pkl'.format(idx))
# files.sort()
run_folder = 'stage3_all/' #!!!
# run_folder = 'stage3_second_half/' #!!!
z_str = str(z)

# filters = ['F150W', 'F356W']
filters = ['F150W']
import copy, matplotlib
for top_psf_id in range(3):
# for top_psf_id in [0]:
    fit_run_list = []
    # idx = idx_info
    filt = filters[0]
    
    if filt == 'F150W' :
        cmap = 'inferno'
    else:
        cmap = 'gist_heat'
    my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap)) # copy the default cmap
    my_cmap.set_bad('black')

    PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
    # idx, filt= item
    # fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
    fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
    # fit_files = [fit_files[i] for i in range(len(fit_files)) if '_IncludeOtherID' not in fit_files[i] and 'useidx' not in fit_files[i]]
    fit_files.sort()
    for i in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort()  
    print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
    fit_run = fit_run_list[sort_Chisq[top_psf_id]]
    # fit_run.savename = 'figures/' + fit_run.savename+'_'+filt
    fit_run.plot_final_qso_fit(target_ID = target_id+'$-$'+filt, save_plot = False, cmap = my_cmap)
    count_n = 5
    Chisq_best = chisqs[sort_Chisq[top_psf_id]]
    Chisq_last= chisqs[sort_Chisq[count_n-1]]
    inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
    weight = np.zeros(len(chisqs))
    for i in sort_Chisq[:count_n]:
        weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
    
    
    prop_name = 'magnitude'
    all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
    # all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
    weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
    rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
    
    
#     # result.append([filt, fit_run.fitting_specify_class.zp, weighted_value, rms_value])
    host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
    AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
    ratio = host_flux/(host_flux+AGN_flux)
    print("HOST ratio", ratio)
    print(fit_run.final_result_galaxy[0]['R_sersic'])
    print(fit_run.final_result_galaxy)
    print(prop_name, round(weighted_value,2), '+-', round(rms_value,2))
#     print('Chisqs top 2', round(chisqs[sort_Chisq[0]],2), round(chisqs[sort_Chisq[1]],2))
#     print_s =filt +' ratio: ' + str(round(ratio,2)) + "\n\n\n"
#     print(fit_files[sort_Chisq[0]])
#     print(print_s)
#     # hold = input(print_s)
# # hold = input("idx {0} above, OK?\n\n".format(idx))
    PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
    PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
    
    PSF_RA_DEC_list = np.array(PSF_RA_DEC_list)
    # print( 'The number to the east:',
    #     np.sum(PSF_RA_DEC_list[sort_Chisq[top_psf_id]][0] - PSF_RA_DEC_list[:,0] > 0))
        
        #%%
# noise_map  = fit_run.fitting_specify_class.data_process_class.noise_map
# host = fit_run.flux_2d_out['data-Point Source']
# mask = fit_run.fitting_specify_class.data_process_class.target_mask
# from galight.tools.astro_tools import plt_fits
# plt.imshow(abs(host)/noise_map * mask, origin='lower', vmin = 1, vmax = 10)
# plt.colorbar()
# plt.show()     