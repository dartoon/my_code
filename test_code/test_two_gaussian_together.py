#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 17:10:51 2019

@author: Dartoon

A test of Two gaussian
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt


sig1 = 3
dis1 = np.random.normal(0, sig1, size=1000)

sig2 = 4
dis2 = np.random.normal(0, sig2, size=1000)

dis_mix = dis1+dis2

dis3 = np.random.normal(0, np.sqrt(sig1**2+sig2**2), size=1000)


plt.figure(figsize=(8,6))
plt.hist(dis_mix,histtype=u'step',normed=True,
         label=('dis_mix sample'), linewidth = 2, color='orange')
plt.hist(dis3, histtype=u'step',normed=True,
         label=('dis3 sample'), linewidth = 2, color='green')
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.yticks([])
plt.show()

from scipy import stats
print(stats.ks_2samp(dis_mix, dis3).pvalue)