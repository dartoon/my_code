#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 13:23:13 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from functions_for_result import esti_smass_fixAv
import glob

# folder = '20221103' #metallicity considered with nebular lines
# folder = '20221104' #No lines; mel -1
# folder = '20221105' #No lines; mel -2.5
# folder = '20221105' #No lines; mel -2.5
# folder = '20221106' #Increase AV 1 for obj 1
# folder = '20221109' #fix mel as -1, age as 0.2
# folder = '20221115' #fix mel as -1, age as 0.5, dust = 1.0
folder = '20221120'
rerun = False

mag_result = {'F356W': 23.23, 'F150W': 25.63}

z = 6.40
# idx = 101
idx = 101  #For id 0 
target_id = 'J2236'
band_as_upper = []

Av_list = [0.1,0.3,0.5,0.7,1]
Mel_list = [-1, -0.7, -0.3]

for i, Mel in enumerate(Mel_list):
    for j,Av in enumerate(Av_list):
        if j == 0:
            flag = 0
        else:
            flag = 1
        esti_smass_fixAv(ID = folder+str(idx)+str(i)+str(j), mags_dict = mag_result, z = z, flag = flag, 
                    if_run_gsf=True, band_as_upper = band_as_upper, metallicity=Mel, Av=Av,
                    mag_err=[0.1]*len(mag_result), just_run = False)

# print(idx)
# steller_file = glob.glob('esti_smass/'+folder+str(idx)+'/SFH_*.fits')[0]
# hdul = pyfits.open(steller_file)
# info = hdul[0].header 

# # print(name_list[idx])
# print( "[{0:.1f}$-${1:.1f}]".format(float(info['Mstel_16']), float(info['Mstel_84'])) )
# print( "[{0:.1f}$-${1:.1f}]".format(10**float(info['SFR_16']), 10**float(info['SFR_84'])) )
# print( "[{0:.1f}$-${1:.1f}]".format(10**float(info['T_MW_16']), 10**float(info['T_MW_84'])) )
    
