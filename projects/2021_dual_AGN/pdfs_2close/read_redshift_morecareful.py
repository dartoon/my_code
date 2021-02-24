#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 12:47:42 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
# DR144.4_short.asc
# 011110.89+001709.2

f = open("DR144.4_short.asc","r")
string = f.read()
lines = string.split('\n')   # Split in to \n

def common_member(a, b):     
    a_set = set(a) 
    b_set = set(b) 
    # check length  
    if len(a_set.intersection(b_set)) > 0: 
        return list(a_set.intersection(b_set))
    else: 
        return [] 

def read_z(ID):
    line_0 = [lines[i] for i in range(len(lines)) if ID.split('.')[0] in lines[i].split('  ')[0]]
    line_1 = [lines[i] for i in range(len(lines)) if ID.split('.')[1][3:] in lines[i].split('  ')[0]]
    line = common_member(line_0, line_1)
    if line != [] and len(line) == 1:
        z = float(line[0].split(' ')[-1])
    elif line != [] and len(line) > 1:
        print("Mulit-lines for ID", ID)
        return -99
    else:
        z = -99
    return z

#%%
files = glob.glob('*_HSC-I_proof-2close.pdf')
z_list = []
ID_list = []
for i in range(len(files)):
    file = files[i]
    ID = file.split('_HSC')[0]
    z = read_z(ID)
    if z >0:
        ID_list.append(ID)
        z_list.append(z)