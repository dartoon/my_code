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
filters = ['F356W']
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
        count_n = 5
        Chisq_best = chisqs[sort_Chisq[top_psf_id]]
        Chisq_last= chisqs[sort_Chisq[count_n-1]]
        inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
        weight = np.zeros(len(chisqs))
        for i in sort_Chisq[:count_n]:
            weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
        
#%% Test fix n as 1-4 or 0.5-2.5
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import copy
# # fit_sepc = fit_run.fitting_specify_class
# _data_process = copy.deepcopy(fit_run.fitting_specify_class.data_process_class)
# fit_sepc = FittingSpecify(_data_process)
# fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 5)
#                           # ps_pix_center_list = [ps_pos]  ) #, fix_n_list= [[0,4],[1,1]])
# fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
# fit_sepc.kwargs_params['lens_light_model'][3][0]['n_sersic'] = 0.5
# fit_sepc.kwargs_params['lens_light_model'][4][0]['n_sersic'] = 2.5
# fit_sepc.plot_fitting_sets()
# fit_run = FittingProcess(fit_sepc, savename = target_id)
# fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
# fit_run.plot_final_qso_fit(target_ID =target_id)
#%% Test fix n as 1-4 or 0.5-2.5
data_process = copy.deepcopy(fit_run.fitting_specify_class.data_process_class)
data_process.target_stamp = fit_run.flux_2d_out['data-point source']
fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 0, supersampling_factor = 5)
                          # ps_pix_center_list = [ps_pos]  ) #, fix_n_list= [[0,4],[1,1]])
fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
# fit_sepc.kwargs_params['lens_light_model'][3][0]['n_sersic'] = 1
# fit_sepc.kwargs_params['lens_light_model'][4][0]['n_sersic'] = 1
fit_sepc.plot_fitting_sets()
fit_run = FittingProcess(fit_sepc, savename = target_id)
fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
fit_run.plot_final_galaxy_fit(target_ID =target_id)
