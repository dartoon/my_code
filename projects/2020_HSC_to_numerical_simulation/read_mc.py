#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 10:14:47 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob
files = glob.glob('offset_result/Hor*')
for i in range(len(files)):
    print(i, files[i])
    
read_i = int(input())
file = files[read_i]

data = np.loadtxt(file)

if file[-7:-4] == '1.5':
    data[:,0] = -data[:,0]
    data[:,1] = -data[:,1]
    data[:,2] = -data[:,2]

cals0 = data[:,0][abs(data[:,0])!=99]
cals1 = data[:,1][abs(data[:,1])!=99]
cals2 = data[:,2][abs(data[:,2])!=99]                  
    
print(file)
print('zs:', file[-7:-4])
print("noise no sl {0:.2f}$\pm${1:.2f}".format(np.mean(cals0), np.std(cals0)) )
print("sim {0:.2f}$\pm${1:.2f}".format(np.mean(cals1), np.std(cals1)) )
print("obs {0:.2f}$\pm${1:.2f}".format(np.mean(cals2), np.std(cals2)) )