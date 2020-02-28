#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 10:10:48 2019

@author: Dartoon

Test the burning in of the MCMC
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

import corner
mat.rcParams['agg.path.chunksize'] = 10000

# folder = 'MCsteps1000_10000_1.3.0_non-linear_True_2ndrun/'
folder = 'final_confirmed_MCsteps1000_10000/'
# folder = 'MCsteps1000_10000_1.3.0_non-linear_False/'

model = 'SPEMD'
#model = 'SIE'

for i in range(1,2):
    result = pickle.load(open(folder+'sampler_results_{1}#{0}.pkl'.format(i, model),'rb'),encoding="latin1") 
    fit_result, trans_result = result
    mcmc_new_list, labels_new = trans_result

kwargs_result, chain_list_mcmc, chain_list_pso = fit_result
sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list_mcmc[0]
chain_num = len(mcmc_new_list)
#%%
print(len(mcmc_new_list[1]),labels_new)

idx = [i for i in range(len(labels_new)) if '$\\phi_{lens}$' == labels_new[i]][0]
del labels_new[idx]
mcmc_new_list = np.delete(mcmc_new_list, idx, axis=1)

if model == 'SIE':
	truths=[1.23651, 0.9, 3898, 70.7]
elif model == 'SPEMD':
	truths=[1.23651, 2.0, 0.9, 3898, 70.7]	

plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True, #range= [[0.8,1.5],[1,3],[0,1],[0, 1],[2000,5000],[20,100]], 
                     quantiles=[0.16, 0.5, 0.84], smooth = 0.4,truths =truths,
                     title_kwargs={"fontsize": 15}, label_kwargs = {"fontsize": 25},
                     levels=1.0 - np.exp(-0.5 * np.array([1.,2.]) ** 2))
plt.savefig("corner_plot_{0}_ID1.pdf".format(model))
plt.show()
