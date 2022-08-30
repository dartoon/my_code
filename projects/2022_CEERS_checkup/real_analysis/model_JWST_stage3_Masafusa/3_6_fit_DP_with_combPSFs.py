#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 17:03:39 2022

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
import sys

dp_files = glob.glob('fit_material/data_process_idx1_F200W_CombPsfsNO_*.pkl')
dp_files.sort()
f = open("../model_JWST_stage3_Masafusa/target_idx_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n

# i = int(sys.argv[1]) - 1 
for i in range(len(dp_files)):
    file = dp_files[i]
    print(file)
    idx = file.split('idx')[1].split('_')[0]
    target_id = [lines[i].split(' ')[1] for i in range(len(lines)) if lines[i].split(' ')[0] == str(idx)][0]
    _data_process = pickle.load(open(file,'rb'))
    _data_process.noise_map = np.nan_to_num(_data_process.noise_map, nan=1000)
    if int(idx) in [31, 32, 56]:
        ps_pos = np.where(_data_process.target_stamp == _data_process.target_stamp.max())
        ps_pos = (ps_pos[0][0] - _data_process.radius, ps_pos[1][0] - _data_process.radius)
        ps_pos = [ps_pos[1], ps_pos[0]]
    else:
        ps_pos = _data_process.apertures[0].positions - _data_process.radius
    
    filename = dp_files[i].replace('data_process', 'fit_run')[:-4]+'_{0}.pkl'.format(i)
    fit_sepc = FittingSpecify(_data_process)
    fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3,
                                  ps_pix_center_list = [ps_pos]  ) #, fix_n_list= [[0,4],[1,1]])
    fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
    fit_sepc.build_fitting_seq()
    # fit_sepc.plot_fitting_sets()
    fit_run = FittingProcess(fit_sepc, savename = target_id, fitting_level=['norm','deep','deep'])
    fit_run.run(algorithm_list = ['PSO','PSO','PSO'])
    fit_run.plot_final_qso_fit(target_ID =target_id)
    # filt = _data_process.filt
    pickle.dump(fit_run , open(filename, 'wb'))