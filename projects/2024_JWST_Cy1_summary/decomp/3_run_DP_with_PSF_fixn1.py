#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 19:30:47 2024

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle, glob

from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import sys
sys.path.insert(0, '../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info


filt = 'F356W'
point_source_supersampling_factor = 2 #!!!
# info = target_info[str(idx)]
# target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
take_folder = '../material/{0}/'.format(filt) #!!!
save_folder = '../material/fit_result/'
supersampling_factor = 3

# idx = 0 #!!!
for idx in range(10):
# for idx in [1]:
    if idx == 1:
        continue
    
    info = target_info[str(idx)]
    target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
    data_process = pickle.load(open(take_folder+'data_process_idx{0}.pkl'.format(idx),'rb'))
    PSF_library_DP_files = glob.glob(take_folder+"PSF_library*_DPformat.*pkl")
    PSF_library_DP_files.sort()
    
    for PSFs_file in PSF_library_DP_files:
        
        psf_idx = PSFs_file.split('idx')[1].split('_')[0]
        psf_list_dp = pickle.load(open(PSFs_file,'rb'))
        PSF_list_clean = psf_list_dp.PSF_list_clean
        for i in range(len(PSF_list_clean)):
            if_run = glob.glob(save_folder+'fit_run_{0}_idx{1}_psfidx{2}_{3}_psfsf{4}_fixn1.pkl'.format(filt, idx, psf_idx, i, point_source_supersampling_factor))
            print(idx, psf_idx, i, point_source_supersampling_factor)
            if if_run != []:
                continue
            psf = PSF_list_clean[i]
            psf[psf<0] = 0.
            psf = abs(psf)
            data_process.PSF_list = [psf]
            data_process.noise_map = np.nan_to_num(data_process.noise_map, nan=1000)
            fit_sepc = FittingSpecify(data_process)
            fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = supersampling_factor,
                                         point_source_supersampling_factor = point_source_supersampling_factor,
                                         apertures_center_focus=True)
            fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
            fit_sepc.kwargs_params['lens_light_model'][4][0]['R_sersic'] = 1.
            fit_sepc.kwargs_params['lens_light_model'][2][0]['n_sersic'] = 1.
            fit_sepc.plot_fitting_sets()
            fit_run = FittingProcess(fit_sepc, savename = target_id)
            fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
            fit_run.plot_final_qso_fit(target_ID =target_id)
            pickle.dump(fit_run , open(save_folder+'fit_run_{0}_idx{1}_psfidx{2}_{3}_psfsf{4}_fixn1.pkl'.format(filt, idx, psf_idx, i, point_source_supersampling_factor), 'wb'))
            