#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 18:30:42 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from matplotlib.ticker import AutoMinorLocator
from scipy import stats

zs = 0.5
files = glob.glob("offset_result/*{0}.txt".format(zs))

if zs == 0.5 or zs == 0.7:
    files = files+glob.glob("offset_result/*{0}.txt".format(0.6))
files.sort()

ct = 0
fig, ax = plt.subplots(figsize=(8,7))

colors = ['green','steelblue','c','deeppink','plum','m']
for i in [3,2,1,4,5,0]:
    data = np.loadtxt(files[i]) #nosl, sl, obs
    if zs == 1.5:
        data = -data
    if ct == 0:
        plt.hist(data[:,2][abs(data[:,2])!=99],histtype=u'step',density=True,
                  label=' Observation', linewidth = 4)
        
    plt.hist(data[:,1][abs(data[:,1])!=99], histtype=u'step',density=True,
              label=' '+files[i].split("/")[1].split('_')[0], linewidth = 2,alpha = 0.8)
              # color = colors[i])
    ct += 1
    pvalue = stats.ks_2samp(data[:,2][abs(data[:,2])!=99], data[:,1][abs(data[:,1])!=99] ).pvalue
    if pvalue>1.e-10:
        print(zs, files[i].split("/")[1].split('_')[0], "%e"%pvalue )
    else:
        print(zs, files[i].split("/")[1].split('_')[0], "<1e-10")

plt.title(r"Offset compared with the local, z={0}".format(zs), fontsize = 25)
plt.tick_params(labelsize=20)
# plt.legend(prop={'size':20})
if zs == 1.5:
    plt.legend(scatterpoints=1,numpoints=1,loc=1,prop={'size':20},ncol=1,handletextpad=0)
plt.yticks([])
# plt.xlim([-2.5,3.5])
ax.xaxis.set_minor_locator(AutoMinorLocator())
plt.tick_params(which='both', width=2, top=True,direction='in')
plt.tick_params(which='major', length=10)
plt.tick_params(which='minor', length=6)#, color='r’)
plt.xlabel(r'$\Delta$log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.savefig('offset_dis_z{0}.pdf'.format(zs))
plt.show()
    