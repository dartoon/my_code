#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 21:48:10 2019

@author: Dartoon
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner
import astropy.io.fits as pyfits
from matplotlib.colors import LogNorm
import copy

band_seq = ['G', 'R', 'I', 'Z', 'Y']
# picklename = '084710.40-001302.6/fit_image_084710.40-001302.6_HSC-I.pkl'
picklename = '121405.12+010205.1/fit_image_121405.12+010205.1_HSC-I.pkl'


result = pickle.load(open(picklename,'rb'))
best_fit, chain_list_result, trans_paras, material = result

source_result, image_host, ps_result, image_ps, _ =best_fit
# chain_list, param_list, _ = pso_fit
# samples_mcmc, param_mcmc, dist_mcmc, _ = mcmc_fit
chain_list, _ = chain_list_result
sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[1]

print("best-fit source_result:", source_result)
print("best-fit ps_result:", ps_result)
pix_scale = 0.167
multi_band_list, kwargs_model, kwargs_result, QSO_msk, kwargs_fixed_source, kwargs_fixed_ps, kwargs_constraints, kwargs_numerics = material

#%%diagnose the PSO chain convergency
#import lenstronomy.Plots.output_plots as out_plot
from lenstronomy.Plots import chain_plot
for i in range(len(chain_list)):
    f, axes = chain_plot.plot_chain_list(chain_list,0)

#%%test the MCMC chain convergency
import lenstronomy.Plots.chain_plot as out_plot
#        
#import lenstronomy.Plots.output_plots as plot_mcmc_behaviour
fig = plt.figure(figsize=(20, 15))
ax = fig.add_subplot(111)
out_plot.plot_mcmc_behaviour(ax, samples_mcmc, param_mcmc, dist_mcmc)       

#%% Recover the plot
from lenstronomy.Plots.model_plot import ModelPlot
modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result,
                          arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[QSO_msk])
f, axes = plt.subplots(3, 3, figsize=(16, 16), sharex=False, sharey=False)
modelPlot.data_plot(ax=axes[0,0], text="Data")
modelPlot.model_plot(ax=axes[0,1])
modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)

modelPlot.decomposition_plot(ax=axes[1,0], text='Host galaxy', source_add=True, unconvolved=True)
modelPlot.decomposition_plot(ax=axes[1,1], text='Host galaxy convolved', source_add=True)
modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)

modelPlot.subtract_from_data_plot(ax=axes[2,0], text='Data - Point Source', point_source_add=True)
modelPlot.subtract_from_data_plot(ax=axes[2,1], text='Data - host galaxy', source_add=True)
modelPlot.subtract_from_data_plot(ax=axes[2,2], text='Data - host galaxy - Point Source', source_add=True, point_source_add=True)

f.tight_layout()
plt.show()

# #%%If need to re-plot the corner plot or translate any parameter:
# from lenstronomy.Sampling.parameters import Param
# param = Param(kwargs_model, kwargs_fixed_source=kwargs_fixed_source, kwargs_fixed_ps=kwargs_fixed_ps, **kwargs_constraints)
# kwargs_out = param.args2kwargs(samples_mcmc[i])

#%%Recover the translated cornor plot
QSO_img = multi_band_list[0][0]['image_data']

plt.imshow(QSO_img, origin='low', norm=LogNorm())
for i in range(len(ps_result)):
    obj_x, obj_y = len(QSO_img)/2 - ps_result[i]['ra_image'][0]/pix_scale, len(QSO_img)/2+ps_result[i]['dec_image'][0]/pix_scale
    # print(obj_x, obj_y)
    plt.text(obj_x, obj_y, "QSO{0}".format(i), fontsize=10, color='k')
plt.show()    

plt.imshow(QSO_img, origin='low', norm=LogNorm())
for i in range(len(source_result)):
    obj_x, obj_y = len(QSO_img)/2 - source_result[i]['center_x']/pix_scale, len(QSO_img)/2+source_result[i]['center_y']/pix_scale
    plt.text(obj_x, obj_y, "obj{0}".format(i), fontsize=15, color='k')
plt.show()   
 
mcmc_new_list, labels_new, _ = trans_paras
plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
plt.show()

#%%The host flux for the host flux:
#Read the fitting parameter.
idx = 0
v_l=np.percentile(samples_mcmc[:,idx],16,axis=0)
v_m=np.percentile(samples_mcmc[:,idx],50,axis=0)
v_h=np.percentile(samples_mcmc[:,idx],84,axis=0)
print(param_mcmc[idx], ":", v_l, v_m, v_h)

#For the translated totol flux.
mcmc_new_list = np.asarray(mcmc_new_list)
idx = 2
v_l=np.percentile(mcmc_new_list[:,idx],16,axis=0)
v_m=np.percentile(mcmc_new_list[:,idx],50,axis=0)
v_h=np.percentile(mcmc_new_list[:,idx],84,axis=0)
print(labels_new[idx], ":", v_l, v_m, v_h)
    