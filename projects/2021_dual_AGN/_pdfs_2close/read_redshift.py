#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:53:34 2021

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

def read_z(ID):
    line = [lines[i] for i in range(len(lines)) if ID in lines[i]]
    if line != []:
        z = float(line[0].split(' ')[-1])
    else:
        z = -99
    return z


#%%
# files = glob.glob('*_HSC-I_proof-2close.pdf')

# write_file = open('detected.txt','w') 
# for i in range(len(files)):
#     file = files[i]
#     ID = file.split('_HSC')[0]
#     line = [lines[i] for i in range(len(lines)) if ID in lines[i]]
#     write_file.write("{0}\n".format(line))
# write_file.close()


# z_list = []
# ID_list = []
# for i in range(len(files)):
#     file = files[i]
#     ID = file.split('_HSC')[0]
#     z = read_z(ID)
#     if z >0:
#         ID_list.append(ID)
#         z_list.append(z)  #In total of 272 that have redshfit information.


# #%%
# import shutil
# for i in range(len(z_list)):
#     ID = ID_list[i]
#     file_name =  ID + '_HSC-I_proof-2close.pdf'
#     if z_list[i] <0.2:
#         shutil.move(file_name, 'z_below_0.2/'+file_name)