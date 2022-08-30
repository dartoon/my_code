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
from galight.tools.measure_tools import mask_obj

dp_files = glob.glob('fit_material/data_process_id*F444W*.pkl')
target_id, RA, Dec = 'SDSS_0', 214.8234583, 52.83030556
dp_files.sort()
#%%
# for i in range(33):
# for i in range(33, 66):
# for i in range(66, 99):
# for i in range(99, len(dp_files)):
for i in range(len(dp_files)):
    filename = dp_files[i].replace('data_process', 'fit_run')[:-4]+'_{0}.pkl'.format(i)
    print(filename)
    if glob.glob(filename) != []:
        continue
    file = dp_files[i]
    _data_process = pickle.load(open(file,'rb'))
    psf = _data_process.PSF_list[-1]
    psf[psf<0] = 0.
    psf = abs(psf)
    _data_process.PSF_list[-1] = psf
    _data_process.noise_map = np.nan_to_num(_data_process.noise_map, nan=1000)

    target_stamp = _data_process.target_stamp 
    mask = mask_obj(target_stamp, _data_process.apertures[:1], if_plot=False)
    ps_pos = np.where(target_stamp == np.max(target_stamp * (1-mask[0])))
    ps_pos = (ps_pos[0][0] - _data_process.radius, ps_pos[1][0] - _data_process.radius)
    ps_pos = [ps_pos[1], ps_pos[0]]
    
    fit_sepc = FittingSpecify(_data_process)
    fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3,
                                  ps_pix_center_list = [ps_pos]  ) #, fix_n_list= [[0,4],[1,1]])
    fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
    # fit_sepc.kwargs_constraints['linear_solver'] = False
    fit_sepc.plot_fitting_sets()
    fit_run = FittingProcess(fit_sepc, savename = target_id)
    fit_run.run(algorithm_list = ['PSO','PSO','PSO'], fitting_level=['norm','deep', 'deep'])
    fit_run.plot_final_qso_fit(target_ID =target_id)
    filt = _data_process.filt
    pickle.dump(fit_run , open(filename, 'wb'))
    
    # fit_run.fitting_specify_class.
    host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
    AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
    ratio = host_flux/(host_flux+AGN_flux)