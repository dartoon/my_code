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

run_folder = 'stage3_all/' #!!!
filt = 'F150W'
# dp_files = glob.glob('fit_material_second_smFOV/data_process_idx0_*F150W*.pkl')
# dp_files = glob.glob('fit_material_second_smFOV/data_process_idx0_F356W_FOVpsf2.pkl')
idx = 4
dp_files = glob.glob(run_folder+'fit_material/data_process_idx{1}_{0}_FOVpsf*.pkl'.format(filt,idx))
dp_files.sort()


import sys
sys.path.insert(0, '../model_z6_data_id0/')
from target_info import target_info

info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

#%%
# for i in [0]:
for i in range(len(dp_files)):
# for i in range(8,16):
# for i in range(16,24):
# for i in range(24, 30):
# for i in range(30, len(dp_files)):
    file = dp_files[i]
    print(file)
    filename = dp_files[i].replace('data_process', 'fit_run')[:-4]+'_{0}.pkl'.format(i)
    idx = file.split('idx')[1].split('_')[0]
    target_id = target_id
    _data_process = pickle.load(open(file,'rb'))
    psf = _data_process.PSF_list[-1]
    psf[psf<0] = 0.
    psf = abs(psf)
    _data_process.PSF_list[-1] = psf
    _data_process.noise_map = np.nan_to_num(_data_process.noise_map, nan=1000)
    
    # if int(idx) in [31, 32, 56]:
    #     ps_pos = np.where(_data_process.target_stamp == _data_process.target_stamp.max())
    #     ps_pos = (ps_pos[0][0] - _data_process.radius, ps_pos[1][0] - _data_process.radius)
    #     ps_pos = [ps_pos[1], ps_pos[0]]
    # else:
    #     ps_pos = _data_process.apertures[0].positions - _data_process.radius
    
    fit_sepc = FittingSpecify(_data_process)
    fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3, apertures_center_focus=True)
                                  # ps_pix_center_list = [ps_pos]  ) #, fix_n_list= [[0,4],[1,1]])
    fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
    fit_sepc.kwargs_params['lens_light_model'][4][0]['R_sersic'] = 1.
    # fit_sepc.kwargs_constraints['linear_solver'] = False
    fit_sepc.plot_fitting_sets()
    fit_run = FittingProcess(fit_sepc, savename = target_id)
    fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
    fit_run.plot_final_qso_fit(target_ID =target_id)
    filt = _data_process.filt
    pickle.dump(fit_run , open(filename, 'wb'))
    
    # # fit_run.fitting_specify_class.
    # host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
    # AGN_flux = fit_run.final_result_ps[0]['flux_within_frame'] 
    # ratio = host_flux/(host_flux+AGN_flux)
    
    