#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 13:23:13 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from functions_for_result import esti_smass
import glob

# folder = '20221103' #metallicity considered with nebular lines
# folder = '20221104' #No lines; mel -1
# folder = '20221105' #No lines; mel -2.5
# folder = '20221105' #No lines; mel -2.5
# folder = '20221106' #Increase AV 1 for obj 1
# folder = '20221109' #fix mel as -1, age as 0.2
# folder = '20221115' #fix mel as -1, age as 0.5, dust = 1.0
# folder = '20221231' #fix mel as -1, age as 0.5, dust = 1.0
folder = '20230119' #fix mel as -1, age as 0.5, dust = 1.0
rerun = False

mag_result = {'F150W': 26.4, 'F356W': 24.80} #QSO host

# mag_result = {'F150W': 22.85, 'F356W': 21.49} #obj1
# mag_result = {'F150W': 24.64, 'F356W': 23.55} #obj2
# mag_result = {'F150W': 26.23, 'F356W': 26.36} #obj3
# del mag_result['F410M']

idx = 2
target_id = 'J2255'
band_as_upper = ['F150W']

# idx, target_id = 201,  'obj1'
# idx, target_id = 202,  'obj2'
# idx, target_id = 203,  'obj3'
z = 6.34
import time
t1 = time.time()
print('Run estimate')
esti_smass(ID = folder+str(idx), mags_dict = mag_result, z = z, flag = 1, 
            if_run_gsf=True, band_as_upper = band_as_upper, metallicity=-0.7,
            mag_err=[0.15]*len(mag_result), just_run = False)
t2 = time.time()

print(idx)
steller_file = glob.glob('esti_smass/'+folder+str(idx)+'/SFH_*.fits')[0]
hdul = pyfits.open(steller_file)
info = hdul[0].header 

# print(name_list[idx])
print( "[{0:.1f}$-${1:.1f}]".format(float(info['Mstel_16']), float(info['Mstel_84'])) )
print( "[{0:.1f}$-${1:.1f}]".format(10**float(info['SFR_16']), 10**float(info['SFR_84'])) )
print( "[{0:.1f}$-${1:.1f}]".format(10**float(info['T_MW_16']), 10**float(info['T_MW_84'])) )
    
    # print('redshift', float(info['ZMC_50']))
    # print('smass',  float(info['Mstel_16']) , float(info['Mstel_50']) , float(info['Mstel_84']) )
    # print('sfr', round(10**float(info['SFR_50']),3) )
    # print('age', round(10**float(info['T_MW_50']),3) )
    # print('AV', float(info['AV_50']) )
    # print('SSFR', round(float(info['SFR_50']) - float(info['Mstel_50']) ,3) )
    # print("open",'esti_smass/'+folder+str(idx), '/' )
    # print('\n')


# # 
# steller_file = glob.glob('esti_smass/'+folder+str(idx)+'/gsf_spec_*.fits')[0]
# hdul = pyfits.open(steller_file)
# info = hdul[0].header