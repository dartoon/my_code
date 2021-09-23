#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 23:36:39 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

CIV = 1549
MgII = 2798
Hb = 4861
Ha = 6563
OIII = 5007

#Cal HST filter
G102range = [8000, 11500]
G141range = [10750, 17000]
z = 1.4444
def av_filter(z):
    lines = np.array([CIV, MgII, Hb, OIII])
    redshift_lines = (1+z) * lines
    G102_bool =  (redshift_lines>G102range[0]+100) * (redshift_lines<G102range[1]-100)
    G141_bool =  (redshift_lines>G141range[0]+100) * (redshift_lines<G141range[1]-100)
    # return G102_bool, G141_bool
    s1 = np.array(['CIV', 'MgII', 'Hb', 'OIII'])[G102_bool] 
    s2 = np.array(['CIV', 'MgII', 'Hb', 'OIII'])[G141_bool] 
    s1 = [s1[i] for i in range(len(s1))]
    s2 = [s2[i] for i in range(len(s2))]
    # str1 = "G102: " + repr(s1)
    # str2 = " G141: " + repr(s2)    
    if s2 == []:
        s = "G102 &" + repr(s1)
    elif s1 != []:
        s = "G141 &" + repr(s2)
    else:
        s = "No fileter!!! & "
    return s

print(av_filter(1.4444))