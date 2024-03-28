#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 13:23:13 2022

@author: Dartoon

Note: The generated files are moved to '../../../2024_JWST_Cy2/'
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
# folder = '20221107' #Increase AV 1 for obj 1
# folder = '20221108' #FIX n =1 result
# folder = '20221120' #FIX n =1 result
# folder = '20230113' #FIX n =1 result, test mel = 0.2
# folder = '20230527' #FIX n =1 result, update F150 result
folder = '20240126' #FIX n =1 result, update Cy2 result
folder = '20240311' #FIX n =4 result, update Cy2 result
# idx = [1,2,0,51,35]

# mag_result = {'F356W': 23.12, 'F150W': 25.12}
# mag_result = {'F356W': 23.12, 'F150W': 25.12, 'F444W': 23.21, 'F480M': 23.23, 'F300M': 23.41, 'F250M': 24.56, 'F200W': 24.68, 'F115W': 25.7}
mag_result = {'F356W': 22.572, 'F150W': 24.902, 'F115W': 25.266, 'F250M': 24.44, 'F444W': 22.628, 'F200W': 24.369, 'F300M': 22.591, 'F480M': 22.725} 
# del mag_result['F410M']

target_id = 'J2236'
# idx = 101 #no emission line, fix mel as -1
# idx = 101 #no emission line, fix mel as -1
# idx = 102 #With emission line, fix mel as -1
# idx = 103 #no emission line, fix mel as -2.5
# idx = 104 #With emission line, fix mel as -2.5
# idx, target_id = 201,  'obj1'
# idx, target_id = 202,  'obj2'
# idx, target_id = 203,  'obj3'
idx = 1 #no emission line, fix mel as -1

z = 6.4
import time
t1 = time.time()
print('Run estimate')
band_as_upper = []
esti_smass(ID = folder+str(idx), mags_dict = mag_result, z = z, flag = 0, 
            if_run_gsf=True, band_as_upper = band_as_upper, metallicity=None,
            mag_err=[0.2, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2], just_run = False)
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