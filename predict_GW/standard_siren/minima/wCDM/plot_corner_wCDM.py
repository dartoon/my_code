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

name = 'wCDM_(5, 20)_5-20'
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
#data=data[data[:,2]<88]
data=data[data[:,0]>0.20]
samples = data
fig = corner.corner(samples, labels=["$\Omega_{m}$","$w$", "$H_0$"],
                    quantiles=[0.16, 0.84],show_titles=True,
                    title_kwargs={"fontsize": 12},truths=[0.3,-1,70],
                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2),
                    title_fmt='.2f', range=[(0.15,0.40),(-1.9,-0.5),(60,90)] )
fig.savefig("wcdm_{0}.pdf".format(name[9:11]))
#####
plt.show()   
'''

################load the data##############
value = 2
import pickle
ndim =3
samplerchain=pickle.load(open("mcmc_wcdm_%s"%(value),'rb'))
burn=samplerchain[:,:,:].T
plt.plot(burn[0,20:,:], '-', color='k', alpha=0.3)  #show the chain after 50 steps 
samples = samplerchain[:, 40:, :].reshape((-1, ndim))
fig = corner.corner(samples, labels=["$\Omega_{m}$","$w$", "$H_0$"],
                    quantiles=[0.16, 0.84],show_titles=True,
                    title_kwargs={"fontsize": 12},truths=[0.3,-1,70],
                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2),
                    title_fmt='.2f', range=[(0.15,0.40),(-1.9,-0.5),(60,90)] )
#fig.savefig("triangle.png")
#####
plt.show()   
'''