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
#filename = 'test2_select-eff_correct_sigmalogdiv3_20.txt'
#filename = 'few_tests/right_answer_test1_conv_lognorm_20.txt'
filename = 'test3_result/test3_mode1_level20.txt'
#print filename
#with open(filename) as f:
#        content = f.readlines()
#lines = [x.strip() for x in content] 
#import re
#lines = [re.findall("\d+\.\d+", lines[i]) for i in range(len(lines))]
#lines = [lines[i] for i in range(len(lines)) if len(lines[i]) ==4]
#numbers = np.asarray(lines).astype(float)
numbers = np.loadtxt(filename)
samples = numbers#[numbers[:,1]<3]
fig = corner.corner(samples, labels=["a0", "a1", "mbh_max", "mbh_min"],
                    truths=[2.35, 0.7, 80., 5.],
                    quantiles=[0.16, 0.84],show_titles=True,
                    title_kwargs={"fontsize": 12},#truths=[2.35,80,5],
#                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2),
                    title_fmt='.2f')
#fig.savefig("test1_lognorm_likeli.pdf")
#####
plt.show()  
