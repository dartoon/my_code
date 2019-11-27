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
#truths=[1.6, 0.2, 50.0, 5.]
#truths=[1.6, 0.7, 50.0, 5.]
#truths=[1.6, 1.2, 50.0, 5.]

truths=[2.4, 0.2, 50.0, 5.]

f = open("201911_newrun/model1_a0{0}_a1{1}_mbhmax-50.0_noizl-20.txt".format(truths[0],truths[1]),"r")
string = f.read()
lines = string.split('\n')   # Split in to \n
result = []

for i in range(len(lines)-1):
    line = lines[i].split('([')[-1]
    line = line.split('])')[0]
    if len(line)>2:
        line = [float(stri) for stri in line.split(',')]
        result.append(line)
        
result = np.asarray(result)

#
#result = np.loadtxt('201911_newrun/model1_a01.6_a10.2_mbhmax-50.0_noizl-20.txt')
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
                    truths=truths,
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
print len(numbers)