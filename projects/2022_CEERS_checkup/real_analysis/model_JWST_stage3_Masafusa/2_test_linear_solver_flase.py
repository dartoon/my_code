#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 14:39:50 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
file = 'fit_material/fit_run_idx10_F115W_CombPsfsNO_8_10.pkl'
fit_run = pickle.load(open(file,'rb'))
# fit_run.plot_final_qso_fit()
print(fit_run.final_result_ps[0]['point_amp'])
print(fit_run.fitting_seq.likelihoodModule.log_likelihood(verbose=True, kwargs_return=fit_run.kwargs_result))

fit_run.fitting_specify_class.kwargs_params['lens_light_model'][0][0]['amp'] = 1
fit_run.fitting_specify_class.kwargs_params['lens_light_model'][0][1]['amp'] = 1
fit_run.fitting_specify_class.kwargs_params['lens_light_model'][1][0]['amp'] = 1
fit_run.fitting_specify_class.kwargs_params['lens_light_model'][1][1]['amp'] = 1
fit_run.fitting_specify_class.kwargs_params['lens_light_model'][3][0]['amp'] = 0
fit_run.fitting_specify_class.kwargs_params['lens_light_model'][3][1]['amp'] = 0
fit_run.fitting_specify_class.kwargs_params['lens_light_model'][4][0]['amp'] = 1.e8
fit_run.fitting_specify_class.kwargs_params['lens_light_model'][4][1]['amp'] = 1.e8

fit_run.fitting_specify_class.kwargs_params['point_source_model'][0][0]['point_amp'] = [1]
fit_run.fitting_specify_class.kwargs_params['point_source_model'][1][0]['point_amp'] = [1]
fit_run.fitting_specify_class.kwargs_params['point_source_model'][3][0]['point_amp'] = [0]
fit_run.fitting_specify_class.kwargs_params['point_source_model'][4][0]['point_amp'] = [1.e8]

fit_sepc = fit_run.fitting_specify_class
from lenstronomy.Workflow.fitting_sequence import FittingSequence
fit_run.fitting_specify_class.fitting_seq = FittingSequence(fit_sepc.kwargs_data_joint, 
                                                            fit_sepc.kwargs_model, 
                                                            fit_sepc.kwargs_constraints, fit_sepc.kwargs_likelihood, 
                                                            fit_sepc.kwargs_params, mpi=fit_sepc.mpi)
fit_run.fitting_specify_class.kwargs_constraints['linear_solver'] = False

fit_run.fitting_level = 'norm'
fit_run.run(algorithm_list = ['PSO','MCMC'])
fit_run.plot_final_qso_fit()

# #%%
# from lenstronomy.Plots.model_plot import ModelPlot
# modelPlot = ModelPlot(fit_run.fitting_specify_class.kwargs_data_joint['multi_band_list'],
#                       fit_run.fitting_specify_class.kwargs_model, fit_run.kwargs_result,
#                       arrow_size=0.02, cmap_string="gist_heat", 
#                       image_likelihood_mask_list=fit_run.fitting_specify_class.kwargs_likelihood['image_likelihood_mask_list'],
#                       )    
    
# f, axes = plt.subplots(3, 3, figsize=(16, 16), sharex=False, sharey=False)
# modelPlot.data_plot(ax=axes[0,0], text="Data")
# modelPlot.model_plot(ax=axes[0,1])
# modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)

# modelPlot.decomposition_plot(ax=axes[1,0], text='Host galaxy', lens_light_add=True, unconvolved=True)
# modelPlot.decomposition_plot(ax=axes[1,1], text='Host galaxy convolved', lens_light_add=True)
# modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)

# modelPlot.subtract_from_data_plot(ax=axes[2,0], text='Data - Point Source', point_source_add=True)
# modelPlot.subtract_from_data_plot(ax=axes[2,1], text='Data - host galaxy', lens_light_add=True)
# modelPlot.subtract_from_data_plot(ax=axes[2,2], text='Data - host galaxy - Point Source', lens_light_add=True, point_source_add=True)
# f.tight_layout()
# plt.show()
