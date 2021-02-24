#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 22:34:17 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import os
import shutil

#Move file to ID_folders
# files = glob.glob('*-I.fits')
# for i in range(len(files)):
#     file = files[i]
#     ID = file.split('_HSC')[0]
#     # os.mkdir('{0}'.format(ID))
#     move_files = glob.glob('*{0}*'.format(ID))
#     for move_file in move_files:
#         shutil.move(move_file, ID)

#%%
f = open("../2close/DR144.4_short.asc","r")
string = f.read()
lines = string.split('\n')   # Split in to \n

def read_z(ID):
    line = [lines[i] for i in range(len(lines)) if ID in lines[i]]
    if line != []:
        z = float(line[0].split(' ')[-1])
    else:
        z = -99
    return z

IDs = glob.glob('??????.*.*')
for ID in IDs:
    if read_z(ID)>1:
         shutil.move(ID, 'z_over1')

