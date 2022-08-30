#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 16:47:18 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob
import pickle
import copy

from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess

import shutil
from galight.tools.measure_tools import mask_obj
#%%
# filt = 'F444W'
# files = glob.glob('fit_material/fit_run_idx1_{0}*.pkl'.format(filt))
# files.sort()
# file = files[0]
# fit_run = pickle.load(open(file,'rb') )


data_process = fit_run.fitting_specify_class.data_process_class
aperture = copy.deepcopy(data_process.apertures[0])
data_process.apertures.append(aperture)

target_stamp = data_process.target_stamp 
mask = mask_obj(target_stamp, data_process.apertures[:1], if_plot=False)
ps_pos = np.where(target_stamp == np.max(target_stamp * (1-mask[0])))
ps_pos = (ps_pos[0][0] - data_process.radius, ps_pos[1][0] - data_process.radius)
ps_pos = [ps_pos[1], ps_pos[0]]

fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3,
                              ps_pix_center_list = [ps_pos]  ) #, fix_n_list= [[0,4],[1,1]])
fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.2
# fit_sepc.kwargs_constraints['linear_solver'] = False
fit_sepc.plot_fitting_sets()
fit_run = FittingProcess(fit_sepc)
fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
fit_run.plot_final_qso_fit()
