#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 20:06:03 2022

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

from functions_for_result import load_prop, load_info
from galight.tools.measure_tools import measure_FWHM

import shutil
from galight.tools.measure_tools import mask_obj
idx = 51
fit_run_dict = load_prop(idx, root_folder = './*', prop_name='fit_run')
Reff_diff = []
refit = False
if refit == True:
    for filt in fit_run_dict.keys():
        fit_run = fit_run_dict[filt]
        data_process = fit_run.fitting_specify_class.data_process_class
        mask_aperture = copy.deepcopy(data_process.apertures[0])
        mask_aperture.positions = np.array([data_process.radius, data_process.radius])
        PSF_FWHM = np.mean(measure_FWHM(data_process.PSF_list[0]))
        mask_aperture.a = PSF_FWHM*1.2
        mask_aperture.b = PSF_FWHM*1.2
            
        mask = mask_obj(data_process.target_stamp, [mask_aperture])[0]
        data_process.target_mask = mask * data_process.target_mask
        data_process.apertures = data_process.apertures[:1] 
        data_process.target_stamp = fit_run.flux_2d_out['data-Point Source'] - np.sum(fit_run.image_host_list[1:],axis=0 )
        data_process.plot_materials()
        fit_sepc = FittingSpecify(data_process)
        fit_sepc.prepare_fitting_seq(point_source_num = 0, supersampling_factor = 3 ) #, fix_n_list= [[0,4],[1,1]])
        fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
        # fit_sepc.kwargs_constraints['linear_solver'] = False
        fit_sepc.plot_fitting_sets()
        fit_run_withmask = FittingProcess(fit_sepc)
        fit_run_withmask.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
        fit_run_withmask.plot_final_galaxy_fit()
        print(fit_run.final_result_galaxy[0]['R_sersic'], fit_run_withmask.final_result_galaxy[0]['R_sersic'])
        Reff_diff.append([filt, fit_run.final_result_galaxy[0]['R_sersic'], 
                          fit_run_withmask.final_result_galaxy[0]['R_sersic']])
        pickle.dump(fit_run_withmask, open('refit_host_mask/idx'+str(idx)+filt+'_host_fit.pkl', 'wb'))


#%%
# for idx in [0]:
for idx in [0, 1, 2, 35, 51]:
    target_id, z = load_info(idx)
    fit_run_dict = load_prop(idx, root_folder = './*', prop_name='fit_run')
    dirct_Reff = []
    refit_Reff = []
    prop = 'R_sersic'
    # prop = 'n_sersic'
    # prop = 'magnitude'
    # prop = 'flux_within_frame'
    plt.figure(figsize=(6,6))
    _record = []
    for filt in fit_run_dict.keys():
        refit_file = 'refit_host_mask/idx'+str(idx)+filt+'_host_fit.pkl'
        fit_run = fit_run_dict[filt]
        refit_run = pickle.load(open(refit_file,'rb'))
        plt.scatter(fit_run.final_result_galaxy[0][prop], refit_run.final_result_galaxy[0][prop],
                    label = filt)
        _record.append([fit_run.final_result_galaxy[0][prop], refit_run.final_result_galaxy[0][prop],])
    xmin = np.min(_record) / 1.2
    xmax = np.max(_record) * 1.2
    plt.xlabel(r"direct inferred "+ prop,fontsize=25)
    plt.ylabel(r"refit withmask "+ prop,fontsize=25)
    plt.plot(np.linspace(0.03,xmax), np.linspace(0.03,xmax))
    plt.xlim([xmin,xmax])
    plt.ylim([xmin,xmax])
    plt.tick_params(labelsize=20)
    plt.title(str(idx)+' '+target_id+" z="+str(z), fontsize=17)
    plt.legend(prop={'size':10})
    plt.show()
    