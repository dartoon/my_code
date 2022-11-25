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
import sys
sys.path.insert(0,'../model_z6_data_id0/')

idx = 1
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# files = glob.glob('./*fit_material*sm*/data_process_idx{0}_*_psf*.pkl'.format(idx))
# files.sort()
# run_folder = 'stage3_*/' #!!!
run_folder = 'stage3_all/' #!!!
z_str = str(z)

# filters = ['F150W', 'F356W']
filters = ['F150W']
import copy, matplotlib
for top_psf_id in [0]:
    for count in range(len(filters)):
        fit_run_list = []
        # idx = idx_info
        filt = filters[count]
        
        if filt == 'F150W' :
            cmap = 'inferno'
        else:
            cmap = 'gist_heat'
        my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap)) # copy the default cmap
        my_cmap.set_bad('black')

        PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
        # idx, filt= item
        fit_files = glob.glob(run_folder+'*fit_material*/fit_run_fixn1__idx{0}_{1}_FOV*.pkl'.format(idx, filt))#+\
        # fit_files = glob.glob(run_folder+'*fit_material*/fit_run_withcentralMask_idx{0}_{1}_FOV*.pkl'.format(idx, filt))#+\
        # fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]
        fit_run.savename = 'figures/' + fit_run.savename+'_'+filt
        fit_run.plot_final_qso_fit(target_ID = target_id+'$-$'+filt, save_plot = True, cmap = my_cmap)
