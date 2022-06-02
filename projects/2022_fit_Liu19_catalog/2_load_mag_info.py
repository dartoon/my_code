#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 13:24:22 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from fun_smass_estimate import esti_smass

f_Liu_t1 = open("Liu_files/table1.dat","r")
t1_string = f_Liu_t1.read()
t1_lines = t1_string.split('\n')

f_Liu_t2 = open("Liu_files/table2.dat","r")
t2_string = f_Liu_t2.read()
t2_lines = t2_string.split('\n')

f_Liu_t1 = open("Liu_files/table1.dat","r")
t1_string = f_Liu_t1.read()
t1_lines = t1_string.split('\n')

f = open("table_mag_frameflux.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
Mbh_list = []
z_list = []
# for line in lines[1:-1]:
for line in lines[1:2]:
    results = line.split()
    ID, magG, magR, magI, magZ, magY = results
    info_l = t1_lines[int(ID)-1]
    info = info_l.split(' ')
    info = [s for s in info if s != '']
    ID_load0, SDSS_ID, z = info[0], info[1], info[4]
    
    info2_l = t2_lines[int(ID)-1]
    info2 = info2_l.split(' ')
    info2 = [s for s in info2 if s != '']
    ID_load1, Mbh = info2[0], info2[-2]
    Mbh_list.append(float(Mbh))
    z_list.append(float(z))
    if ID!= ID_load0 or ID!= ID_load1:
        print(ID)
    mags = magG, magR, magI, magZ, magY
    mags = [float(mag) for mag in mags]
    z = float(z)
    esti_smass(ID, mags, z)