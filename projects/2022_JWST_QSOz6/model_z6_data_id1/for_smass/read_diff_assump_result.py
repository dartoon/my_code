#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 23:48:59 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

import sys
sys.path.insert(0,'../')

folder = '20221120' #

fitidx = 0
from target_info import target_info
info = target_info[str(fitidx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

import matplotlib as mat

idx = 101  #Age free
# idx = 102  #Age fix 0.5

folder = 'esti_smass/'+folder+str(idx)
steller_files = glob.glob(folder+'*/gsf_spec_*.fits')

result = []
for steller_file in steller_files:
    hdul = pyfits.open(steller_file)
    info_muv = hdul[1].header 
    # print(info_muv['MUV50'], info_muv['MUV16'], info_muv['MUV84'])
    
    steller_file = glob.glob(steller_file.replace('gsf_spec_', 'SFH_'))[0]
    hdul = pyfits.open(steller_file)
    info1 = hdul[0].header 
    # print('redshift', float(info1['ZMC_50']))
    # print('smass', info1['Mstel_50'], info1['Mstel_16'], info1['Mstel_84']) 
    
    #%%
    mat.rcParams['font.family'] = 'STIXGeneral'
    hdul = pyfits.open(steller_file.replace('gsf_spec_', 'SFH_'))
    info = hdul[0].header 
    # print(target_id)
    sfr = round(10**float(info['SFR_50']),3) 
    z = info['Z']
    sfr_l = round(10**float(info['SFR_16']),3) 
    sfr_h = round(10**float(info['SFR_84']),3) 
    smass_l = round(float(info['Mstel_16']),3) 
    smass = round(float(info['Mstel_50']),3) 
    smass_h = round(float(info['Mstel_84']),3) 
    age = round(10**float(info['T_MW_50']),3) 
    age_l = round(10**float(info['T_MW_16']),3) 
    age_h = round(10**float(info['T_MW_84']),3) 
    mel = round(float(info['Z_MW_50']),1) 
    Av = round(float(info['AV_50']),1) 
    Muv = round(info_muv['MUV50'],2)
    
    # print('Av, mel, Age, M*', info['AV_50'], mel, age, smass)
    result.append([Av, mel, age, smass, Muv])


# Av_list = [0.1,0.3,0.5,0.7,1]
Av_list = [0.3,0.5,0.7]
Mel_list = [-1, -0.7, -0.3]

print('Av, mel, Age, M*, Muv')
for i, Mel in enumerate(Mel_list):
    for j,Av in enumerate(Av_list):
        for k in range(len(result)):
            if result[k][0] == Av and result[k][1] == Mel:
                print(result[k])
                

# print('redshift', float(info['ZMC_50']))
# print('smass', float(info['Mstel_50']) )
# print('sfr', sfr)
# print('sfr_l', sfr_l)
# print('sfr_h', sfr_h)
# print('age', age)
# print('age_l', age_l)
# print('age_h', age_h)
# print('AV', float(info['AV_50']) )
# print('SSFR', round(float(info['SFR_50']) - float(info['Mstel_50']) ,3) )
# # print("open",'esti_smass/'+folder+str(idx), '/' )
# print('\n')

