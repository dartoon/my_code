#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 16:55:50 2024

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

idx = 6

info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
fit_run_list = []
filt = 'F150W'
count_n = 5
all_library = True

# for add_cond in ['', '_fixn1']:
for add_cond in ['_fixn1']:
    fit_run_list = []
    if all_library == True:
        load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx*_psfsf*{2}.pkl'.format(filt, idx, add_cond))
    else:
        psf_idx = idx
        load_files = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx{2}*_psfsf*{3}.pkl'.format(filt, idx, psf_idx,add_cond))
    load_files.sort()
    chisqs_idx = []
    for file in load_files:
        fit_run_list.append(pickle.load(open(file,'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort()  
    file_name = load_files[sort_Chisq[0]]   
    fit_run = fit_run_list[sort_Chisq[0]]    
    fit_run.plot_final_qso_fit()


# #%%
# from galight.fitting_specify import FittingSpecify
# from galight.fitting_process import FittingProcess

# rerun = True
# point_source_supersampling_factor  = int(file_name.split('psfsf')[1][0])
# if rerun == True:
#     supersampling_factor = 3
#     fit_sepc = FittingSpecify(fit_run.fitting_specify_class.data_process_class)
#     fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = supersampling_factor,
#                                  point_source_supersampling_factor = point_source_supersampling_factor,
#                                  apertures_center_focus=True)
#     fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
#     fit_sepc.kwargs_params['lens_light_model'][4][0]['R_sersic'] = 10*fit_run.fitting_specify_class.deltaPix
#     fit_sepc.plot_fitting_sets()
#     fit_run_newrun = FittingProcess(fit_sepc, savename = target_id)
#     fit_run_newrun.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
#     fit_run_newrun.plot_final_qso_fit(target_ID =target_id)


