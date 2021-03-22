#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 18:21:50 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

f = open("../_pdfs_2close/DR144.4_short.asc","r")
# f = open("material/ID_RA_DEC_z.txt","r")
string = f.read()
zlines = string.split('\n')   # Split in to \n
# write_file = open('ID_RA_DEC_z.txt','w') 

def read_z(ID):
    line = [zlines[i] for i in range(len(zlines)) if ID in zlines[i]]
    if line != []:
        line  = line[0].split(' ')
        line = [x for x in line if x!='']        
        z = float(line[-1])
    else:
        z = -99
    return z

def read_info(ID):
    line = [zlines[i] for i in range(len(zlines)) if ID in zlines[i]]
    if line != []:
        line  = line[0].split(' ')
        line = [x for x in line if x!='']
        # write_file.write(line[0] + '\n')
        RA = float(line[1])
        Dec = float(line[2])
        z = float(line[-1])
    else:
        z = -99
        RA = - 99
        Dec = - 99
    return RA, Dec, z