#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 16:51:04 2017

@author: Dartoon
"""

'''
See error of the data
'''
import numpy as np
import matplotlib.pyplot as plt

f=open('TL_sim_NO_10000.txt')
data=np.loadtxt(f)
f.close
plt.errorbar(data[:,0],np.log10(data[:,2]),yerr=0.00001,fmt='.',color='blue')#,markersize=1)
plt.xlabel("$z$")
plt.ylabel("Luminosity distance in Mpc, logspace")
plt.show()
