#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 16:22:04 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob


f_Liu_t1 = open("Liu_files/table1.dat","r")
t1_string = f_Liu_t1.read()
t1_lines = t1_string.split('\n')

f_Liu_t2 = open("Liu_files/table2.dat","r")
t2_string = f_Liu_t2.read()
t2_lines = t2_string.split('\n')

f_Liu_t3 = open("Liu_files/table3.dat","r")
t3_string = f_Liu_t3.read()
t3_lines = t3_string.split('\n')

f = open("table_mag_frameflux.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
Mbh_list = []
z_list = []


# write_file = open('table_summary.txt','w') 
# write_file.write("#ID, z, smass, Mbh, Gmag, Rmag, Imag, Zmag, Ymag (based on flux in frame) \n")

# # for line in lines[1:-1]:
# for line in lines[1:-1]:
#     results = line.split()
#     ID, magG, magR, magI, magZ, magY = results
#     info_l = t1_lines[int(ID)-1]
#     info = info_l.split(' ')
#     info = [s for s in info if s != '']
#     ID_load0, SDSS_ID, z = info[0], info[1], info[4]
    
#     info2_l = t2_lines[int(ID)-1]
#     info2 = info2_l.split(' ')
#     info2 = [s for s in info2 if s != '']
#     ID_load1, Mbh = info2[0], info2[-4]
#     print(Mbh)
#     z_list.append(float(z))
#     if ID!= ID_load0 or ID!= ID_load1:
#         print(ID)
#     mags = magG, magR, magI, magZ, magY
#     mags = [float(mag) for mag in mags]
#     z = float(z)
    
#     filename_p  = 'esti_smass/'+ID+'/SFH_*.fits'
#     if glob.glob(filename_p) !=[]:
#         filename_p = glob.glob(filename_p)[0]
#         hdul = pyfits.open(filename_p)
#         table = hdul[1].data
#         name = hdul[1].columns
#         smass_idx = [i for i in range(len(name)) if 'Mstel50' in str(name[i])][0]
#         inf_smass = table[1][smass_idx] # 'smass:'
#         win = [ID, str(z), '{0:.3f} {1:.3f}'.format(inf_smass, float(Mbh)), magG, magR, magI, magZ, magY]
#     else:
#         inf_smass = -99
#         win = [ID, str(z), '{0} {1:.3f}'.format(inf_smass, float(Mbh)), magG, magR, magI, magZ, magY]
#     write = ''
#     for w in win:
#         write = write + w + ' '
#     win = win[0] + win[0]
#     write_file.write(write)
#     # print(ID, z, '{0:.3f} {1:.3f}'.format(inf_smass, float(Mbh)), magG, magR, magI, magZ, magY)
#     write_file.write('\n')
# write_file.close()