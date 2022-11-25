#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:02:58 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
from galight.data_process import DataProcess
from galight.tools.astro_tools import plt_fits
from galight.tools.measure_tools import measure_bkg
import pickle

import pickle

idx = 0
# filters = ['F150W', 'F356W']
filters = ['F356W']
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# files = glob.glob('./*fit_material*sm*/data_process_idx{0}_*_psf*.pkl'.format(idx))
# files.sort()
run_folder = '../model_z6_data_id{0}/stage3_all/'.format(idx) #!!!
z_str = str(z)

import copy, matplotlib
for top_psf_id in [0]:
    for count in range(len(filters)):
        fit_run_list = []
        # idx = idx_info
        filt = filters[count]

        PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
        if idx ==0:
            fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        elif idx ==1:
            fit_files = glob.glob(run_folder+'*fit_material*/fit_run*_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        # print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]

#%%
_data_process = fit_run.fitting_specify_class.data_process_class
psf = _data_process.PSF_list[-1]
psf[psf<0] = 0.
psf = abs(psf)
_data_process.PSF_list[-1] = psf
_data_process.noise_map = np.nan_to_num(_data_process.noise_map, nan=1000)

_data_process.apertures = []
fit_sepc = FittingSpecify(_data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3)
fit_sepc.plot_fitting_sets()
fit_run = FittingProcess(fit_sepc, savename = target_id)
fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
fit_run.plot_final_qso_fit(target_ID =target_id)
