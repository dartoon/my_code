#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 22:58:18 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

H0_SIE = [[68.236, 4.597],
[70.779, 5.179],
[153.034, 4.940],
[186.416, 6.478],
[86.195, 1.292],
[261.970, 9.697],
[76.069, 3.316],
[84.799, 4.623],
[158.870, 4.653],
[112.556, 4.025],
[84.071, 1.611],
[172.391, 4.506],
[84.818, 4.278],
[78.842, 4.571],
[147.805, 4.054],
[169.789, 5.808],
[78.086, 1.398],
[226.757, 7.352]]
#H0_SIE = np.array(H0_SIE)

H0_SPEMD = [[86.191, 4.139],
[89.046, 4.172],
[150.843, 6.505],
[168.679, 7.399],
[78.865, 3.457],
[163.620, 6.077],
[63.572, 3.519],
[39.075, 1.378],
[134.131, 5.052],
[142.343, 4.589],
[112.455, 1.897],
[172.769, 4.449],
[67.301, 2.992],
[57.190, 2.293],
[146.317, 6.295],
[137.776, 6.762],
[116.748, 1.591],
[137.854, 4.185]]
#H0_SPEMD = np.array(H0_SPEMD)

for seed in range(0,6):
    H0_s = np.array([H0_SIE[i][0] for i in range(len(H0_SIE)) if i%6 == seed])
    H0_error = np.array([H0_SIE[i][1] for i in range(len(H0_SIE)) if i%6 == seed])
    ave = np.sum(H0_s/H0_error**2)/np.sum(1/H0_error**2)
    std = np.sqrt( np.sum((H0_s-ave)**2/H0_error**2) / np.sum(1/H0_error**2) )
#    print '#',seed+1, "H0=", round(ave,3),'stdd' , round(np.std(H0_s),3)
    print '#',seed+1, "H0=", round(ave,3),'stdd' , round(std,3)
        
print '\n\n'
for seed in range(0,6):
    H0_s = np.array([H0_SPEMD[i][0] for i in range(len(H0_SPEMD)) if i%6 == seed])
    H0_error = np.array([H0_SPEMD[i][1] for i in range(len(H0_SPEMD)) if i%6 == seed])
    ave = np.sum(H0_s/H0_error**2)/np.sum(1/H0_error**2)
    std = np.sqrt( np.sum((H0_s-ave)**2/H0_error**2) / np.sum(1/H0_error**2) )
#    print '#',seed+1, "H0=", round(np.sum(H0_s/H0_error**2)/np.sum(1/H0_error**2),3),'stdd', round(np.std(H0_s),3)
    print '#',seed+1, "H0=", round(ave,3),'stdd' , round(std,3)
        