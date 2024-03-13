#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 16:58:35 2024

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


# Populate the subplots with content and remove ticks
for jj in range(10):
    idx = jj
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
        
    #%%Save host with WCS
    file_header = copy.deepcopy(fit_run.fitting_specify_class.data_process_class.header)
    file_header['CRPIX1'] = file_header['CRPIX1']-fit_run.fitting_specify_class.data_process_class.target_pos[0]+len(fit_run.image_host_list[0])/2
    file_header['CRPIX2'] = file_header['CRPIX2']-fit_run.fitting_specify_class.data_process_class.target_pos[1]+len(fit_run.image_host_list[0])/2
    pyfits.PrimaryHDU(fit_run.flux_2d_out['data-point source'],header=file_header).writeto('data_minus_QSO/{0}_host.fits'.format(target_id),overwrite=True)
    # data_minus_QSO
    