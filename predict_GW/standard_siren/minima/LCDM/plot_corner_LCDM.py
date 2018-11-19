#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 23:23:23 2018

@author: dartoon
"""
import corner
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams["font.size"] = 16

#level=input("which level?: 5? 10? 15? 20?:\n")
#f0=open('minima_LDCM_sc3_{0}'.format(level))

#prior = input('Which sernaio?')
'''
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
                    title_kwargs={"fontsize": 20},truths=[0.3,70],label_kwargs= {"fontsize": 20},
                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2),
                    title_fmt='.2f', range=[(0.15,0.45),(60,90)] )
fig.savefig("lcdm_{0}.pdf".format(name[9:11]))
####
plt.show()   
'''

################load the data##############
value = 10
import pickle
ndim =2
samplerchain=pickle.load(open("mcmc_lcdm_%s"%(value),'rb'))
burn=samplerchain[:,:,:].T
plt.plot(burn[0,20:,:], '-', color='k', alpha=0.3)  #show the chain after 50 steps 
samples = samplerchain[:, 40:, :].reshape((-1, ndim))
import corner
fig = corner.corner(samples, labels=["$\Omega_{m}$", "$H_0$"],
                    quantiles=[0.16, 0.84],show_titles=True,
                    title_kwargs={"fontsize": 12},truths=[0.3,70],
                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2),
                    title_fmt='.2f', range=[(0.15,0.45),(60,90)] )
#fig.savefig("triangle.png")
#####
plt.show()
h0_l = np.percentile(samples, 16,axis=0)[1]
h0_h = np.percentile(samples, 84,axis=0)[1]
h0_m = np.percentile(samples, 50,axis=0)[1]
level = (h0_h - h0_l)/h0_m/2
print "H0 level:", level*100
