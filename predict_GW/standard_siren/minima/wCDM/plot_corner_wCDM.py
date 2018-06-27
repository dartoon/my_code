#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 23:23:23 2018

@author: dartoon
"""
import corner
import numpy as np
import matplotlib.pyplot as plt

level=input("which level?: 5? 10? 15? 20?:\n")
f0=open('minima_wDCM_{0}_1st'.format(level))
data=np.loadtxt(f0)
data=data[data[:,0]>0]
orl=len(data)
#data=data[data[:,1]>-1.98]
#data=data[data[:,2]<94]
#data=data[data[:,0]<0.54]
samples = data
fig = corner.corner(samples, labels=["$\Omega_{m}$","$w$", "$H_0$"],
                    quantiles=[0.16, 0.5, 0.84], title_kwargs={"fontsize": 12},truths=[0.3,-1,70],
                    plot_datapoints=False,smooth=1.0,smooth1d=1.0,plot_density=True, levels=(0.6826, 0.9544),
                    title_fmt='.3f',show_titles=False )
#fig.savefig("wcdm_%s.pdf"%(level))
#####
plt.show()   

