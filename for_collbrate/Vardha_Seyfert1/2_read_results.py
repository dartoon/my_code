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

picklename = 'dump_fitting.pkl'
result = pickle.load(open(picklename,'rb'))
[source_result, image_host, ps_result, image_ps, samples_mcmc, param_mcmc, paras] = result

print "best-fit source_result:", source_result
print "best-fit ps_result:", ps_result
#%%

# here the (non-converged) MCMC chain of the non-linear parameters
if not samples_mcmc == []:
    plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
    plt.show()
    [source_params_2, ps_param_2, mcmc_new_list, labels_new] = paras
    plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
    plt.show()

#%% To readout the parameter

#Read the fitting parameter.
idx = 0
v_l=np.percentile(samples_mcmc[:,idx],16,axis=0)
v_m=np.percentile(samples_mcmc[:,idx],50,axis=0)
v_h=np.percentile(samples_mcmc[:,idx],84,axis=0)
print param_mcmc[idx], ":", v_l, v_m, v_h

#For the translated totol flux.
mcmc_new_list = np.asarray(mcmc_new_list)
idx = 2
v_l=np.percentile(mcmc_new_list[:,idx],16,axis=0)
v_m=np.percentile(mcmc_new_list[:,idx],50,axis=0)
v_h=np.percentile(mcmc_new_list[:,idx],84,axis=0)
print labels_new[idx], ":", v_l, v_m, v_h