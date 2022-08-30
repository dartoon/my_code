#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 19:29:07 2022

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

files = glob.glob('fit_material/data_process_idx*.pkl')
files.sort()
collect_info = []
i_list = []
for i in range(len(files)):
# for i in [2768]:
    _file = files[i]
    idx_info = _file.split('idx')[1].split('_')[0]
    filt_info = _file.split('_psf')[0].split('_')[-1]
    if int(idx_info) in [31, 32, 56]:
        i_list.append(i)

# dp_files = glob.glob('fit_material/data_process_idx*.pkl')
# dp_files.sort()
# f = open("target_idx_info.txt","r")
# string = f.read()
# lines = string.split('\n')   # Split in to \
# filt = ''
# import sys
# count = int(sys.argv[1]) - 1

# for i in i_list[count:count+1]:
#     file = dp_files[i]
#     # print(file)
#     filename = dp_files[i].replace('data_process', 'fit_run')[:-4]+'_{0}.pkl'.format(i)
#     idx = file.split('idx')[1].split('_')[0]
#     target_id = [lines[i].split(' ')[1] for i in range(len(lines)) if lines[i].split(' ')[0] == str(idx)][0]
#     _data_process = pickle.load(open(file,'rb'))
#     # psf = _data_process.PSF_list[0]
#     # psf[psf<0] = 0
#     # _data_process.PSF_list[0] = psf
#     _data_process.noise_map = np.nan_to_num(_data_process.noise_map, nan=1000)
#     if filt !=  _data_process.filt:
#         filt = _data_process.filt
#         # ps_pos = _data_process.apertures[0].positions - _data_process.radius
#         ps_pos = np.where(_data_process.target_stamp == _data_process.target_stamp.max())
#         ps_pos = (ps_pos[0][0] - _data_process.radius, ps_pos[1][0] - _data_process.radius)
#         fit_sepc = FittingSpecify(_data_process)
#         fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3,
#                                       ps_pix_center_list = [[ps_pos[1], ps_pos[0]]]  ) #, fix_n_list= [[0,4],[1,1]])
#         fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
#         fit_sepc.build_fitting_seq()
#         print(idx, filt)
#         # fit_sepc.plot_fitting_sets()
#     fit_run = FittingProcess(fit_sepc, savename = target_id, fitting_level=['norm','deep','deep'])
#     fit_run.run(algorithm_list = ['PSO','PSO','PSO'])
#     # fit_run.plot_final_qso_fit(target_ID =target_id)
#     pickle.dump(fit_run , open(filename, 'wb'))
    
#     # # fit_run.fitting_specify_class.
#     # host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
#     # AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
#     # ratio = host_flux/(host_flux+AGN_flux)