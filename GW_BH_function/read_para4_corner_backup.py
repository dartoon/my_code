#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 13:52:19 2019

@author: Dartoon

Read corner
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import matplotlib
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
cmap = matplotlib.cm.get_cmap('viridis')

import corner
#filename = 'test2_select-eff_correct_sigmalogdiv3_20.txt'
#filename = 'few_tests/right_answer_test1_conv_lognorm_20.txt'
result = np.loadtxt('test3_result_secondrun/test3_mode1_take2_level20.txt')
numbers = result


#%%
#print filename
#with open(filename) as f:
#        content = f.readlines()
#lines = [x.strip() for x in content] 
#import re
#lines = [re.findall("\d+\.\d+", lines[i]) for i in range(len(lines))]
#lines = [lines[i] for i in range(len(lines)) if len(lines[i]) ==4]
#numbers = np.asarray(lines).astype(float)
samples = numbers#[numbers[:,1]<3]
fig = corner.corner(numbers, labels=[r"$\alpha_0$", r"$\alpha_1$", "M$_{max}$", "M$_{min}$"],
                    truths=[2.35, 0.7, 80., 5.],
                    quantiles=[0.16, 0.5, 0.84],show_titles=True, smooth = 0.7,
                    title_kwargs={"fontsize": 15}, label_kwargs = {"fontsize": 25},
#                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
                    levels=1.0 - np.exp(-0.5 * np.array([1.,2.]) ** 2),
                    title_fmt='.2f')
for ax in fig.get_axes():
      ax.tick_params(axis='both', labelsize=20)
#fig.savefig("fig_results_4para.pdf")
#####
plt.show()  
