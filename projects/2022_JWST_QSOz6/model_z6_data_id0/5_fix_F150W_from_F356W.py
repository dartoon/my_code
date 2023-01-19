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

import sys
sys.path.insert(0,'../model_z6_data_id0')

run_folder = 'stage3_all/' #!!!

idx = 0
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

from galight.tools.measure_tools import mask_obj
import copy, matplotlib

mag_list = []
chisq_list = []
fit_run_list_fix = []
for top_psf_id in range(5):
    for filt in ['F356W', 'F150W']:
        fit_run_list = []
        fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_FOV*.pkl'.format(idx, filt))#+\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort() 
        if filt == 'F356W':
            fit_run_F356W = fit_run_list[sort_Chisq[0]]
        elif filt == 'F150W':
            fit_run_F150W = fit_run_list[sort_Chisq[top_psf_id]]

    fit_sepc = copy.deepcopy(fit_run_F150W.fitting_specify_class)
    fit_sepc.kwargs_params['lens_light_model'][2][0]['R_sersic'] = fit_run_F356W.final_result_galaxy[0]['R_sersic']
    fit_sepc.kwargs_params['lens_light_model'][1][0]['R_sersic'] = fit_run_F356W.final_result_galaxy[0]['R_sersic']
    fit_sepc.kwargs_params['lens_light_model'][2][0]['n_sersic'] = fit_run_F356W.final_result_galaxy[0]['n_sersic']
    fit_sepc.kwargs_params['lens_light_model'][1][0]['n_sersic'] = fit_run_F356W.final_result_galaxy[0]['n_sersic']
    fit_sepc.kwargs_params['lens_light_model'][2][0]['e1'] = fit_run_F356W.final_result_galaxy[0]['e1']
    fit_sepc.kwargs_params['lens_light_model'][1][0]['e1'] = fit_run_F356W.final_result_galaxy[0]['e1']
    fit_sepc.kwargs_params['lens_light_model'][2][0]['e2'] = fit_run_F356W.final_result_galaxy[0]['e2']
    fit_sepc.kwargs_params['lens_light_model'][1][0]['e2'] = fit_run_F356W.final_result_galaxy[0]['e2']
    
    fit_run = FittingProcess(fit_sepc, savename = target_id)
    fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
    fit_run.plot_final_qso_fit(target_ID =target_id)
    fit_run_F150W_fix = fit_run
    fit_run_list_fix.append(fit_run_F150W_fix)
    mag_list.append(fit_run_F150W_fix.final_result_galaxy[0]['magnitude'])
    chisq_list.append(fit_run_F150W_fix.reduced_Chisq)

    #%%

chisqs = chisq_list
Chisq_best = np.min(chisq_list)
Chisq_last= np.max(chisq_list)
inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
weight = np.zeros(len(chisqs))
for i in range(len(chisq_list)):
    weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
    
# all_values = mag_list
# weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
# rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))

prop_name = 'magnitude'
# all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
all_values = [fit_run_list_fix[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list_fix))]
weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
print('host', prop_name, "{0:.2f}$\pm${1:.2f}".format(weighted_value, rms_value))

prop_name = 'magnitude'
all_values = [fit_run_list_fix[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list_fix))]
weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
print('quasar', prop_name, "{0:.2f}$\pm${1:.2f}".format(weighted_value, rms_value))

prop_name = 'flux_within_frame'
all_values = [100*(fit_run_list_fix[i].final_result_galaxy[0][prop_name]/ (fit_run_list_fix[i].final_result_galaxy[0][prop_name] + fit_run_list[i].final_result_ps[0][prop_name]))
              for i in range(len(fit_run_list_fix))]
weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
print('host flux ratio " ' , "{0:.1f}\%$\pm${1:.1f}\%".format(weighted_value, rms_value))
