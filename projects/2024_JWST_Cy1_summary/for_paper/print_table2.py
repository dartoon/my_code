#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 14:37:23 2024

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
from galight.tools.plot_tools import plot_data_apertures_point
from galight.tools.astro_tools import plt_fits
from photutils.aperture import EllipticalAperture
import sys
import copy, matplotlib

import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, '../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info

run_folder = '../material/fit_result/'
prop_name = 'magnitude'
# prop_name = 'R_sersic'
# prop_name = 'n_sersic'
count_n = 5 
# add_cond = ''  #!!!

# idx = 0   
filt = ['F356W','F150W'][1] #!!!
for idx in range(10 ):
    if idx == 0:
        all_library = False #!!! for the new run that use all the PSF library, or just only in FOV.
    else:
        all_library = True #!!! for the new run that use all the PSF library, or just only in FOV.
    info = target_info[str(idx)]
    target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
    run_folder_2022 = '../../2022_JWST_QSOz6/model_z6_data_id{0}/stage3_all/'.format(idx) #!!!
    z_str = str(z)
    
    target_id = target_id.replace('-', '$-$')
    print_s = target_id + ' & '
    
    add_cond_list = ['_fixn1', '']
    for count_i, add_cond in enumerate(add_cond_list):
        # add_cond = '_fixn1'  #!!!
        fit_run_list_sp1 = []
        psf_sp = 1
        if all_library == True:
            load_files_sp1 = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx*_psfsf{2}{3}.pkl'.format(filt, idx, psf_sp, add_cond))
        else:
            psf_idx = idx
            load_files_sp1 = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx{2}*_psfsf{3}{4}.pkl'.format(filt, idx, psf_idx, psf_sp,add_cond))
        load_files_sp1.sort()
        chisqs_idx = []
        for file in load_files_sp1:
            fit_run_list_sp1.append(pickle.load(open(file,'rb')))
        chisqs = np.array([fit_run_list_sp1[i].reduced_Chisq for i in range(len(fit_run_list_sp1))])
        sort_Chisq_sp1 = chisqs.argsort()  
        weight_sp1 = np.zeros(len(chisqs))
        for i in sort_Chisq_sp1[:count_n]:
            weight_sp1[i] = 1
        psf_sp = 2
        fit_run_list_sp2 = []
        if all_library == True:
            load_files_sp2 = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx*_psfsf{2}{3}.pkl'.format(filt, idx, psf_sp, add_cond))
        else:
            psf_idx = idx
            load_files_sp2 = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx{2}*_psfsf{3}{4}.pkl'.format(filt, idx, psf_idx, psf_sp,add_cond))
        load_files_sp2.sort()
        chisqs_idx = []
        for file in load_files_sp2:
            fit_run_list_sp2.append(pickle.load(open(file,'rb')))
        chisqs = np.array([fit_run_list_sp2[i].reduced_Chisq for i in range(len(fit_run_list_sp2))])
        sort_Chisq_sp2 = chisqs.argsort()  
        weight_sp2 = np.zeros(len(chisqs))
        for i in sort_Chisq_sp2[:count_n]:
            weight_sp2[i] = 1
        weight = np.concatenate([weight_sp1, weight_sp2])
        fit_run_list = fit_run_list_sp1 + fit_run_list_sp2
        #%%get the host ratio
        if add_cond == '_fixn1':  #!!!:
            fit_run = fit_run_list_sp1[sort_Chisq_sp1[0]]
            fit_run_name = load_files_sp1[sort_Chisq_sp1[0]]
            if fit_run.reduced_Chisq > fit_run_list_sp2[sort_Chisq_sp2[0]].reduced_Chisq:
                fit_run = fit_run_list_sp2[sort_Chisq_sp2[0]]  #Save Top result
                fit_run_name = load_files_sp2[sort_Chisq_sp2[0]]
            all_host_flux = [fit_run_list[i].final_result_galaxy[0]['flux_within_frame'] for i in range(len(fit_run_list))] 
            all_AGN_flux = [fit_run_list[i].final_result_ps[0]['flux_within_frame'] for i in range(len(fit_run_list))] 
            all_values = [all_host_flux[i]/(all_host_flux[i] + all_AGN_flux[i]) for i in range(len(all_host_flux)) ]
            weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
            rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
            print_s = print_s + '{0:.1f}\%$\\pm${1:.1f}\%'.format(weighted_value*100, rms_value*100) + ' & '
        
        #%%get host magnitude
        prop_name = 'magnitude'
        all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))] 
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print_s = print_s + '{0:.2f}$\\pm${1:.2f}'.format(weighted_value, rms_value) + ' & '
        
        #%%get host magnitude
        prop_name = 'R_sersic'
        all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))] 
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print_s = print_s + '{0:.2f}$\\pm${1:.2f}'.format(weighted_value, rms_value) + ' & '
    
        #%%Check Chisq when not add host Sersic
        if count_i == 1:
            psfid1, psfid2 = fit_run_name.split('psfidx')[1].split('_')[:2]
            load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx{2}_{3}_psfsf*_noHost.pkl'.format(filt, idx, psfid1, psfid2 ))
            load_files.sort()
            chisqs_idx = []
            fit_run_list_noHost = []
            for file in load_files:
                fit_run_list_noHost.append(pickle.load(open(file,'rb')))
            chisqs = np.array([fit_run_list_noHost[i].reduced_Chisq for i in range(len(fit_run_list_noHost))])
            best_chisqs  = np.min(chisqs)
            print_s = print_s + '{0:.2f}\%'.format(-(fit_run.reduced_Chisq - best_chisqs)/best_chisqs * 100) + '&'
        
            # if idx == 3:
            #     fit_run.plot_final_qso_fit(target_ID = target_id)
            #     fit_run_check = fit_run
        
    used_PSF_list = [fit_run_list[i].fitting_specify_class.data_process_class.PSF_list[0] for i in range(len(fit_run_list)) if weight[i] != 0]
    PSF_noise_map = np.std(used_PSF_list, axis=0) * fit_run.final_result_ps[0]['flux_within_frame']
    ct = int(abs( len(fit_run.fitting_specify_class.data_process_class.noise_map) - len(PSF_noise_map))/2)
    total_noise = np.sqrt(PSF_noise_map[ct:-ct, ct:-ct] ** 2 +  fit_run.fitting_specify_class.data_process_class.noise_map**2)
    host_signal = fit_run.flux_2d_out['data-point source'] - np.sum(fit_run.image_host_list[1:],axis = 0)

    fit_run.cal_astrometry()
    x = fit_run.final_result_galaxy[0]['position_xy'][0] + len(fit_run.image_host_list[0])/2
    y = fit_run.final_result_galaxy[0]['position_xy'][1] + len(fit_run.image_host_list[0])/2
    a = fit_run.final_result_galaxy[0]['R_sersic']/fit_run.fitting_specify_class.deltaPix * 4
    
    if a > len(fit_run.image_host_list[0])/2 * 0.8:
        a = len(fit_run.image_host_list[0])/2 * 0.8
    
    b = a*fit_run.final_result_galaxy[0]['q']
    theta = - fit_run.final_result_galaxy[0]['phi_G']
    aperture = EllipticalAperture((x, y), a, b, theta=theta)
    SNR_map = host_signal/total_noise
    # if idx == 3:
    #     plt_fits(SNR_map,colorbar=True)
    # plot_data_apertures_point(host_signal , [aperture])
    SNR = np.sqrt(aperture.do_photometry( SNR_map**2 )[0][0])
    print_s = print_s+ '{0:.2f}'.format(SNR) + '\\\\'
    print(print_s)
    aperture = EllipticalAperture((x, y), a, b, theta=theta)
    # print(aperture.do_photometry(host_signal/total_noise**2)[0][0] / aperture.do_photometry(1/total_noise**2)[0][0] / np.sqrt(1/aperture.do_photometry(1/total_noise**2)[0][0]) )
    
