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
name = 'LCDM_(5, 20)_5-20'
f = open(name,"r")
with f as g:
    lines = g.readlines()

data = []
for i in range(1,len(lines)-1):
    line = lines[i].split(' ')
    data.append([float(line[0]),float(line[1])])

data = np.asarray(data)
#data=data[data[:,0]>0]
#data=data[data[:,0]<0.5499]
samples = data
#fig = corner.corner(samples, labels=["$\Omega_{m}$", "$H_0$"],
#                    quantiles=[0.16, 0.5, 0.84],show_titles=True,
#                    title_kwargs={"fontsize": 12},truths=[0.3,70],
#                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,plot_density=True, levels=(1-np.exp(-0.5),1-np.exp(-4*0.5)),
#                    title_fmt='.2f')

fig = corner.corner(samples, labels=["$\Omega_{m}$", "$H_0$"],
                    quantiles=[0.16, 0.84],show_titles=True,
                    title_kwargs={"fontsize": 12},truths=[0.3,70],
                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2),
                    title_fmt='.2f', range=[(0.15,0.45),(60,90)] )
fig.savefig("lcdm_{0}.pdf".format(name[9:12]))
####
plt.show()   

