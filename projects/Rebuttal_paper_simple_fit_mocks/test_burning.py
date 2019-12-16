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

model = 'SPEMD'
for i in range(6,7):
    result = pickle.load(open('Final_before_submitpaper/Final-MCMC_4000_steps/sampler_results_{1}#{0}.pkl'.format(i, model),'rb'))
    fit_result, trans_result = result
    mcmc_new_list, labels_new = trans_result

kwargs_result, chain_list_mcmc, chain_list_pso = fit_result
sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list_mcmc[0]
chain_num = len(mcmc_new_list)

print "plot how value change with chains:"
print labels_new
i = 3
plt.figure(figsize=(11, 8))
plt.plot(range(chain_num), mcmc_new_list[:,i], '-' ,linewidth = 0.1, color='gray')
plt.xlabel("number of Chains",fontsize=27)
plt.ylabel(labels_new[i],fontsize=27)
plt.tick_params(labelsize=20)
plt.show()

print "plot how value change with chains:"
print param_mcmc
i = 2
plt.figure(figsize=(11, 8))
plt.plot(range(len(samples_mcmc)), samples_mcmc[:,i], '-',linewidth = 0.01, color='gray')
plt.xlabel("number of Chains",fontsize=27)
plt.ylabel(param_mcmc[i],fontsize=27)
plt.tick_params(labelsize=20)
plt.show()