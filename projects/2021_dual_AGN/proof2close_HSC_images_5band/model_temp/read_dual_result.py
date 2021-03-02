#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 15:08:15 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
ID = '124618.51-001750.2'
folder = '../z_over1/' + ID + '/fit_result/'
files = glob.glob(folder + 'fit_result_*txt')

band = ['G', 'R', 'I', 'Z', 'Y']
for i in range(len(files)):
    l = [ j for j in range(len(files)) if band[i]+'-band' in files[j]]
    if l != []:
        file = files[l[0]]
        f = open(file,"r")
        string = f.read()
        lines = string.split('\n')   # Split in to \n
        l0 = [ j for j in range(len(lines)) if 'AGN mag:' in lines[j]]
        AGN_mag = lines[l0[1]]
        AGN_mag = AGN_mag.split(' ')[2:4]
        l1 = [ j for j in range(len(lines)) if 'PS PS' in lines[j]]
        offset = lines[l1[1]].split(' ')[-1]
        print('AGN mags', band[i]+'-band', AGN_mag[0], AGN_mag[1], '; pos offset:', offset)
    else:
        print(band[i]+'-band not fitting.')
    