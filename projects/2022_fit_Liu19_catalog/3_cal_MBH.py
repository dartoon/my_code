#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 15:16:12 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

# def cal_MBH(logLHadr, FWHM_a , logL5100dr , FWHM_b):
#     cal_logMa = 6.71+0.48*(logLHadr-42)+2.12*np.log10(FWHM_a/1000)  # as used in Andreas#!!!
#     cal_logMb = 6.91+0.5*(logL5100dr-44)+2.*np.log10(FWHM_b/1000)  # as used in Andreas #!!!
#     return cal_logMa, cal_logMa, cal_logMb

f_Liu_t2 = open("Liu_files/table2.dat","r")
t2_string = f_Liu_t2.read()
t2_lines = t2_string.split('\n')


write_file = open('table_cal_MBH.txt','w') 
write_file.write("#ID, Mha, Mhb \n")

for line in t2_lines[:-1]:
    line = line.split(' ')
    info2 = [s for s in line if s != '']
    ID = info2[0]
    logLHadr = float(info2[3])
    FWHM_a = float(info2[5])
    logL5100dr = float(info2[58])
    FWHM_b = float(info2[37])
    cal_logMa = 6.71+0.48*(logLHadr-42)+2.12*np.log10(FWHM_a/1000)  # as used in Andreas#!!!
    cal_logMb = 6.91+0.5*(logL5100dr-44)+2.*np.log10(FWHM_b/1000)  # as used in Andreas #!!!
    print(cal_logMa-float(info2[-3]), cal_logMb-float(info2[-4]))
    if not cal_logMa>0:
        cal_logMa = -99
    if not cal_logMb>0:
        cal_logMb = -99
    write = ID + ' ' + str(cal_logMa) + ' ' + str(cal_logMb)
    write_file.write(write)
    write_file.write('\n')
write_file.close()
    
    