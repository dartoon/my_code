#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 09:16:35 2018

@author: Dartoon
"""
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
#import sys
#sys.path.insert(0,'../')
from pz import pz
from gene_data import gene_data

import matplotlib as matt
import matplotlib.lines as mlines
from matplotlib import colors
mat.rcParams['font.family'] = 'STIXGeneral'

# =============================================================================
# Test differnet cosmology
# =============================================================================
fig = plt.figure(figsize=(15,6))
ax = fig.add_subplot(111)
err_level= input("Input the error level:\n")
dl = gene_data(10000,err_level, H0=70., om=0.3)
#dl_wr =  gene_data(10000,3, H0=57., om=0.7)
##dl_wwr =  gene_data(10000, 20, H0=65., om=0.54)
ax.errorbar(x=dl[:,0],y=dl[:,2], yerr=dl[:,3],fmt='b.',zorder=-1)
ax.set_yscale('log')

plt.title('The $D_L$ vs redshift, {0} per cent level'.format(err_level),fontsize=15)
plt.xlabel('$z$',fontsize=15)
plt.ylabel('$D_L$', fontsize=15)
plt.tick_params(labelsize=15)
plt.show()