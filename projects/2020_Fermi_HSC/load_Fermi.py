#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 13:41:54 2020

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits


# filename_p  = 'catalog/gll_psc_v20.fit'
filename_p  = 'catalog/gll_psc_v27.fit'
hdul = pyfits.open(filename_p)
table = hdul[1].data
id_name = hdul[1].columns

#bll BL Lac type of blazar
#fsrq FSRQ type of blazar

idx  = [i for i in range(len(id_name)) if 'CLASS1' in str(id_name[i])][0]

# t_type = 'bll'
t_type = 'fsrq'

# filename_ascii = 'Fermi_' + t_type +'.txt'
# _ascii =  open(filename_ascii,'w')
data, pos_list, id_list = [], [], []
for i in range(len(table)):
    if table[i][idx] == t_type:
        data.append(table[i])
        id_list.append(table[i][0])
        # pos_list.append([table[i][2], table[i][3]])  #RA, DEC
        pos_list.append([table[i][-4], table[i][-3]])  #counterpart RA, DEC
#         ### _ascii.write('{0} {1} {2}\n'.format(table[i][0].split(' ')[1], table[i][2], table[i][3]))
#         _ascii.write('{0} {1} {2}\n'.format(table[i][0].split(' ')[1], table[i][-4], table[i][-3])) #counterpart RA, DEC
# _ascii.close()

#%% Generate the exist catalog:
import glob
if t_type == 'bll':
    files = glob.glob('data/BL_Lac/*_HSC-I.fits')
else:
    files = glob.glob('data/FSRQ/*_HSC-I.fits')
filename_ascii = 'Fermi_' + t_type +'_HSC_exist.txt'
_ascii =  open(filename_ascii,'w')
_ascii =  open(filename_ascii,'w')
for i in range(len(files)):
    ID = files[i].split('/')[-1].split('_')[0]
    line = [j for j in range(len(id_list)) if ID in id_list[j]][0]
    _ascii.write('{0} {1} {2}\n'.format(ID, round(pos_list[line][0],6), round(pos_list[line][1],6) ) )
_ascii.close()
