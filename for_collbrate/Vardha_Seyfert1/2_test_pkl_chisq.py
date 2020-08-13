#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 21:48:10 2019

@author: Dartoon
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner
import astropy.io.fits as pyfits
from matplotlib.colors import LogNorm
import copy
import lenstronomy.Util.param_util as param_util

picklename = 'zoutput/' + 'l22_sersicdisk_1.pkl'
picklename = 'zoutput/' + 'l24_sersicdisk_1.pkl'
picklename = 'zoutput/' + 'l35_sersicdisk_1.pkl'

result = pickle.load(open(picklename,'rb'))
best_fit, chain_list_result, trans_paras, material = result
source_result, image_host, ps_result, image_ps, _ =best_fit
chain_list, _ = chain_list_result
if chain_list[-1][0] == 'EMCEE':
    sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[-1]
    mcmc_new_list, labels_new, _ = trans_paras
multi_band_list, kwargs_model, kwargs_result, QSO_msk, kwargs_fixed_source, kwargs_fixed_ps, kwargs_constraints, kwargs_numerics, classes = material

#%% Recover the plot
from lenstronomy.Plots.model_plot import ModelPlot
modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result,
                          arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[QSO_msk])
# f, axes = plt.subplots(3, 3, figsize=(16, 16), sharex=False, sharey=False)
# modelPlot.data_plot(ax=axes[0,0], text="Data")
# modelPlot.model_plot(ax=axes[0,1])
# modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
# modelPlot.decomposition_plot(ax=axes[1,0], text='Host galaxy', source_add=True, unconvolved=True)
# modelPlot.decomposition_plot(ax=axes[1,1], text='Host galaxy convolved', source_add=True)
# modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
# modelPlot.subtract_from_data_plot(ax=axes[2,0], text='Data - Point Source', point_source_add=True)
# modelPlot.subtract_from_data_plot(ax=axes[2,1], text='Data - host galaxy', source_add=True)
# modelPlot.subtract_from_data_plot(ax=axes[2,2], text='Data - host galaxy - Point Source', source_add=True, point_source_add=True)
# f.tight_layout()
# plt.show()

logL = modelPlot._imageModel.likelihood_data_given_model(source_marg=False, linear_prior=None, **kwargs_result)
n_data = modelPlot._imageModel.num_data_evaluate
print(- logL * 2 / n_data, 'reduced X^2 of all evaluated imaging data combined.')

model, _, _, _= modelPlot._imageModel.image_linear_solve(inv_bool=True, **kwargs_result)
model = model[0]
data = multi_band_list[0][0]['image_data']
noise_map = multi_band_list[0][0]['noise_map']
chisq_map = (data - model)**2/noise_map**2 * QSO_msk

plt.imshow( chisq_map ,origin='lower',norm=LogNorm())
plt.colorbar()
plt.show()

print(np.average(chisq_map) )
print('numbers of pixel have chisq larger than 1000:', chisq_map[chisq_map>1000].shape[0] )
chisq_map[chisq_map>1000] = 1 # remove these pixels to 1
print(np.average(chisq_map) )

#%%
from lenstronomy.LightModel.light_model import LightModel
light = LightModel(kwargs_model['source_light_model_list'])
total_flux = light.total_flux(kwargs_result['kwargs_source'])
