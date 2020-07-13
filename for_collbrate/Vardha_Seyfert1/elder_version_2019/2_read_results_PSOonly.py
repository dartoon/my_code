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

picklename = 'result_vardha/PSO_only/l24_sersicdisk.pkl'
picklename = 'result_vardha/l5548_sersicdisk.pkl'
picklename = 'result_vardha/l79_sersicdiskbar.pkl'

result = pickle.load(open(picklename,'rb'))
[best_fit,pso_fit,mcmc_fit, trans_paras] = result

source_result, image_host, ps_result, image_ps, _ =best_fit
chain_list, param_list, _ = pso_fit
samples_mcmc, param_mcmc, dist_mcmc, _ = mcmc_fit

print("best-fit source_result:", source_result)
print("best-fit ps_result:", ps_result)
pix_sz = 0.04

#%%diagnose the PSO chain convergency
#import lenstronomy.Plots.output_plots as out_plot
import lenstronomy.Plots.chain_plot as out_plot
for i in range(len(chain_list)):
    if len(param_list[i]) > 0:
        f, axes = out_plot.plot_chain(chain_list[i], param_list[i])
plt.show()

print('inferred host flux', [image_host[i].sum() for i in range(len(image_host))])