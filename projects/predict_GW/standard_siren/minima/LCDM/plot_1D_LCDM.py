#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 09:00:25 2018

@author: dartoon
"""

import numpy as np
import matplotlib.pyplot as plt
#level=input("which level?: 5? 10? 15? 20?:\n")
#f0=open('minima_LDCM_sc3_{0}'.format(level))

np.set_printoptions(precision=4)
import matplotlib as mat
import matplotlib.lines as mlines
from matplotlib import colors
mat.rcParams['font.family'] = 'STIXGeneral'
plt.figure(figsize=(13,10))

scenario =5 # input('Which sernaio?')
f0=open('minima_LDCM_sc{0}_fixOM'.format(scenario))
data=np.loadtxt(f0)
samples = data

l =np.percentile(samples,16,axis=0)
m =np.percentile(samples,50,axis=0)
h =np.percentile(samples,84,axis=0)
y = np.linspace(0,700,100)
x_true = y*0+ 70
x_l = y*0+ l
x_m = y*0 + m
x_h = y*0 + h
plt.plot(x_true, y, 'red',linewidth=4.0)
#plt.text(70, 500,'true values', fontsize=25)
plt.plot(x_l, y, 'k--',linewidth=4.0)
plt.plot(x_m, y, 'blue',linewidth=4.0)
plt.plot(x_h, y, 'k--',linewidth=4.0)

plt.axis([69,77,0,600])
plt.hist(samples,bins=30)
plt.xlabel('$H_0$',fontsize=35)
plt.ylabel('PDF',fontsize=35)
plt.tick_params(labelsize=25)
plt.yticks([])


#plt.text(m,650,"$%s %s %s$"%(l,m,h))
#fig.savefig("lcdm_%s.pdf"%(level))
#####
   
plt.savefig("1D-hist.pdf")
plt.show()