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

name = 'minima_wDCM_10'
f = open(name,"r")
with f as g:
    lines = g.readlines()


data = []
for i in range(1,len(lines)-1):
    if 'nan' not in lines[i]:
        line = lines[i].split(' ')
        data.append([float(line[0]),float(line[1]),float(line[2])])
data = np.asarray(data)
#data=data[data[:,0]>0]
#orl=len(data)
data=data[data[:,1]>-1.98]
data=data[data[:,1]<-0.55]
data=data[data[:,2]<94]
data=data[data[:,0]<0.54]
samples = data
fig = corner.corner(samples, labels=["$\Omega_{m}$","$w$", "$H_0$"],
                    quantiles=[0.16, 0.5, 0.84], title_kwargs={"fontsize": 12},truths=[0.3,-1,70],
                    plot_datapoints=False,smooth=1.0,smooth1d=1.0,plot_density=True, levels=(0.6826, 0.9544),
                    title_fmt='.3f',show_titles=False )
#fig.savefig("wcdm_%s.pdf"%(level))
#####
plt.show()   

