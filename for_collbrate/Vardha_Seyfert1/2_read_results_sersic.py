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

# ID = '11' #sys.argv[1]

# rf = open("./fitting_smallpsf.dat",'r')
# for line in rf:
#     ob = line.split()[0]
#     if ob == ID:
#         snratio = int(line.split()[1])
#         numberpixels = int(line.split()[2])
#         pixsz = float(line.split()[3])
#         zero = float(line.split()[4])
#         exptime = int(line.split()[5])
#         fr = line.split()[6]
#         kernel_size = int(line.split()[7])
#         psfname = line.split()[8]
#         if len(line.split())>9:
#             mask = np.array(line.split()[9].split(',')).astype('int')
#         else:
#             mask=[]
#         print (snratio, numberpixels, pixsz, zero, exptime, fr, kernel_size, psfname, mask)

# pix_sz = pixsz

# picklename = 'zoutput/l'+ID+'_sersicdiskbar.pkl'
# picklename = 'zoutput/' + 'l11_sersicdiskbar_6.pkl'

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

logL = modelPlot._imageModel.likelihood_data_given_model(source_marg=False, linear_prior=None, **kwargs_result)
n_data = modelPlot._imageModel.num_data_evaluate
print(- logL * 2 / n_data, 'reduced X^2 of all evaluated imaging data combined.')

"""
#%%diagnose the PSO chain convergency
from lenstronomy.Plots import chain_plot
for i in range(len(chain_list)):
    f, axes = chain_plot.plot_chain_list(chain_list,i)
    plt.show()
        
print ('inferred host flux', [image_host[i].sum() for i in range(len(image_host))])
#%%For MCMC
if chain_list[-1][0] == 'EMCEE':
#If you want to see the corner map for all the fitting parameters:
    plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
    plt.show()
    plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
    plt.show()

#%% Checkout the fitting.
print ("best-fit source_result:", source_result)
print ("best-fit ps_result:", ps_result)

phi0, q0 = param_util.ellipticity2phi_q(source_result[0]['e1'], source_result[0]['e2'])
phi1, q1 = param_util.ellipticity2phi_q(source_result[1]['e1'], source_result[1]['e2'])
phi2, q2 = param_util.ellipticity2phi_q(source_result[2]['e1'], source_result[2]['e2'])

# print ("best_fit source_result phi (bulge, disk, bar): {0:.2f}  {1:.2f}  {2:.2f} ".format(phi0, phi1, phi2))
print ("best_fit source_result q (bulge, disk, bar): {0:.4f}  {1:.4f}  {2:.4f} ".format(q0, q1, q2))

for i in range(len(kwargs_result['kwargs_source'])):
    print("Result of No. {0} component:".format(i))
    print("Reff:", round(kwargs_result['kwargs_source'][i]['R_sersic'],2))
    print("Sersic n:", round(kwargs_result['kwargs_source'][i]['n_sersic'],2))

#%%Print the fitted component separately to see what's going on:
for i in range(len(image_host)):
    plt.imshow(image_host[i], norm = LogNorm(), origin='lower')    
    plt.show()

#%%Test if the prior is not work.
sys.path.insert(0,'./fitting_tools/')    
from fit_qso import condition_bulgediskbar
prior_test = condition_bulgediskbar(kwargs_lens = None, kwargs_source=source_result, kwargs_lens_light = None, kwargs_ps = None, kwargs_special = None, kwargs_extinction = None)
if prior_test < 0:
    print("Warning: The Prior not work!")
"""    