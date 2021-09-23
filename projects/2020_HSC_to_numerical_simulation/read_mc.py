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
files = glob.glob('MC_result/Ill*')
for i in range(len(files)):
    print(i, files[i])
    
read_i = int(input())
file = files[read_i]

data = np.loadtxt(file)
print(file)
print("{0:.2f}$\pm${1:.2f}".format(np.mean(data[:,0]), np.std(data[:,0])) )
print("{0:.2f}$\pm${1:.2f}".format(np.mean(data[:,1]), np.std(data[:,1])) )