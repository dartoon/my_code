#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 23:23:23 2018

@author: dartoon
"""
import corner
import numpy as np
import matplotlib.pyplot as plt

#level=input("which level?: 5? 10? 15? 20?:\n")
#f0=open('minima_LDCM_sc3_{0}'.format(level))

#prior = input('Which sernaio?')
f0 = open ('minima_LDCM_dz_020')#.format(prior))
data=np.loadtxt(f0)
orl=len(data)
data=data[data[:,0]>0]
data=data[data[:,0]<0.5499]
samples = data
fig = corner.corner(samples, labels=["$\Omega_{m}$", "$H_0$"],
                    quantiles=[0.16, 0.5, 0.84],
                    title_kwargs={"fontsize": 12},truths=[0.3,70],
                    plot_datapoints=False,smooth=1.0,smooth1d=1.0,plot_density=True, levels=(0.6826, 0.9544),
                    title_fmt='.3f',show_titles=False )
#fig.savefig("lcdm_%s.pdf"%(level))
#####
plt.show()   

