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

# from acstools import acszpt
dp_files = glob.glob('fit_material/data_process_idx*_F606W*.pkl') + glob.glob('fit_material/data_process_idx*_F814W*.pkl') 
dp_files.sort()

#%%
f = open("../model_JWST_stage3_Masafusa/target_idx_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
# for i in range(1):
# for i in range(len(dp_files)):
for i in range(141):
# for i in range(141, 282):
# for i in range(282, 423):
# for i in range(423, len(dp_files)):
    file = dp_files[i]
    print(file)
    idx = file.split('idx')[1].split('_')[0]
    target_id = [lines[i].split(' ')[1] for i in range(len(lines)) if lines[i].split(' ')[0] == str(idx)][0]
    _data_process = pickle.load(open(file,'rb'))
    filename = dp_files[i].replace('data_process', 'fit_run')[:-4]+'_{0}.pkl'.format(i)
    if np.sum(_data_process.target_mask == False) < 10 and glob.glob(filename) != []:
        continue
    elif np.sum(_data_process.target_mask == False) > 10:
        _data_process.noise_map = np.nan_to_num(_data_process.noise_map, nan=1000)
        
    _data_process.target_mask = _data_process.target_stamp != 0
    ps_pos = _data_process.apertures[0].positions - _data_process.radius
    fit_sepc = FittingSpecify(_data_process)
    fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3,
                                  ps_pix_center_list = [ps_pos]  ) #, fix_n_list= [[0,4],[1,1]])
    fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
    fit_sepc.build_fitting_seq()
    fit_sepc.plot_fitting_sets()
    fit_run = FittingProcess(fit_sepc, savename = target_id, fitting_level=['norm','deep','deep'])
    fit_run.run(algorithm_list = ['PSO','PSO','PSO'])
    fit_run.plot_final_qso_fit(target_ID =target_id)
    filt = _data_process.filt
    # pickle.dump(fit_run , open('fit_material/'+'fit_run_idx{2}_{0}_psf{1}.pkl'.format(filt, i, idx), 'wb'))
    pickle.dump(fit_run , open(filename, 'wb'))