#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:57:45 2021

@author: Dartoon
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

#Information is taken from SDSS website:

f = open("SDSS_info/SDSS_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
meaning = lines[0].split(' ')
table_info = [lines[i] for i in range(len(lines)) if lines[i][0] != '#']
info = [table_info[i].split(' ') for i in range(len(table_info))]

from astropy.coordinates import SkyCoord
from astropy import units as u

SDSS_info = []
RA_DEC_arr = []
for i in range(len(info)):
    RA_h, Dec_h = info[i][5], info[i][6]
    pos = SkyCoord('{0} {1}'.format(RA_h, Dec_h), unit=(u.hourangle, u.deg))
    RA, Dec = pos.ra.degree, pos.dec.degree
    SDSS_info.append([RA, Dec, float(info[i][8]), float(info[i][9]), info[i][11]])
    RA_DEC_arr.append([RA, Dec])
    
RA_DEC_arr = np.array(RA_DEC_arr)    
    
f = open("Fermi_bll_HSC_exist.txt","r")
string = f.read()
lines_bll = string.split('\n')   # Split in to \n

f = open("Fermi_fsrq_HSC_exist.txt","r")
string = f.read()
lines_fsrq = string.split('\n')   # Split in to \n
lines = lines_bll + lines_fsrq

#%%
def find_ind(inp_list, inp_ind, find_list):
    diff = np.sqrt(np.sum((find_list-inp_list[inp_ind])**2,axis=1)) * 3600 # Their difference in arcsec
    out_ind = np.where(diff==diff.min())[0][0]
    if diff.min() < 0.15:
        return out_ind, diff.min()
    else:
        return None, diff.min()

def load_SDSS_info(ID):
    # ID = 'J1245.8+0232'
    info = 'Not exist'
    line = [lines[i] for i in range(len(lines)) if ID in lines[i]][0]
    RA, Dec = line.split(' ')[1:] #Load the RA and Dec in Fermi 
    RA, Dec = float(RA), float(Dec)
    diff = np.sqrt(np.sum((np.array([RA, Dec])-RA_DEC_arr)**2,axis=1)) * 3600 # Their difference in arcsec
    if diff.min() < 10:  #Within 1 arcsec
        info = SDSS_info[np.where(diff == diff.min())[0][0]]
    return info
# test = load_SDSS_info('J1222.5+0414')
#%%
reads = lines_bll
reads = lines_fsrq
for i in range(len(reads)):    
    ID = reads[i].split(' ')[0]
    if load_SDSS_info(ID) != 'Not exist':
        print(ID, load_SDSS_info(ID))