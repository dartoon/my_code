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
import corner

import matplotlib
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
cmap = matplotlib.cm.get_cmap('viridis')

#truth = [0.8, 50.0, 5.]
#truth = [1.6, 50.0, 5.]
#truth = [2.4, 50.0, 5.]
alpha, Mmax, Mmin = truth

filename = '201911_newrun/model0_a-{0}_mbhmax-{1}_noizl-20.txt'.format(alpha, Mmax)
#filename = 'test2_select-eff_correct_sigmalogdiv3_20.txt'
#filename = 'few_tests/right_answer_test1_conv_lognorm_20.txt'
#bool_1 = (pre_result[:,1]>2.7)
#bool_2 = (pre_result[:,3]>95)
#bools = np.invert(bool_1+bool_2)
#result_rm_edge = pre_result[bools]
#result_new_edge = np.loadtxt('test3_result/test3_edge.txt')
#numbers = np.concatenate((result_rm_edge,result_new_edge))
#%%
print(filename)
with open(filename) as f:
        content = f.readlines()
lines = [x.strip() for x in content] 
import re
lines = [re.findall("\d+\.\d+", lines[i]) for i in range(len(lines))]
lines = [lines[i] for i in range(len(lines)) if len(lines[i]) ==3]
numbers = np.asarray(lines).astype(float)
#samples = numbers#[numbers[:,1]<3]
fig = corner.corner(numbers, labels=[r"$\alpha$", r"M$_{\rm max}$", r"M$_{\rm min}$"],
                    truths=truth,
                    quantiles=[0.16, 0.5, 0.84],show_titles=True, smooth = 1.0,  #0.7
                    title_kwargs={"fontsize": 15}, label_kwargs = {"fontsize": 25},
#                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
                    levels=1.0 - np.exp(-0.5 * np.array([1.,2.]) ** 2),
                    title_fmt='.2f')
for ax in fig.get_axes():
      ax.tick_params(axis='both', labelsize=20)
# fig.savefig("3para_contour_a0_{0}.pdf".format(alpha))
#####
plt.show()  
