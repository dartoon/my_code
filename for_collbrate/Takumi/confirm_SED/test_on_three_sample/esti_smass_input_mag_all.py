#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 13:23:13 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

files = glob.glob('./photo_share/*st5.csv')
files.sort()

#%%

import pandas as pd
from functions_for_result import esti_smass, load_prop, load_info, name_list
# for i in range(3):
# for i in range(67, len(files)):
for i in range(7, 8):
    idx = '1129{0}'.format(i)
    idx = int(idx)
    # target_id = ['cid_54', 'cid_361', 'cid_452'][int(str(idx)[-1])]
    filename = files[i]
    if i <=67:
        target_id = 'cid_'+filename.split('cid_')[1].split('_')[0]
    else:
        target_id = 'cw-apr-'+filename.split('cw-apr-')[1].split('_')[0]
    z = filename.split('=')[1].split('_')[0]
    
    # info_file =  glob.glob('*{0}*.csv'.format(target_id))[0]
    df = pd.read_csv(filename)
    
    mag_result = {}
    mag_result_err = []
    band_as_upper = []
    for i in range(len(df)):
        if df['photo'][i] > 0:
            mag_result[df['filternames'][i].upper()] = -2.5*np.log10(df['photo'][i])
        elif df['photo'][i] == 0:
            mag_result[df['filternames'][i].upper()] = -2.5*np.log10(df['photo_unc'][i]*3)
            band_as_upper.append(df['filternames'][i].upper())
    #         # mag_result_err[df['filternames'][i].upper()] = abs( -2.5*np.log10(df['photo'][i]+df['photo_unc'][i]) + 2.5*np.log10(df['photo'][i]))
    #         mag_result_err.append(abs( -2.5*np.log10(df['photo'][i]+df['photo_unc'][i]) + 2.5*np.log10(df['photo'][i])))
    
    
    rerun = True
    
    # idx = 101
    target_id = str(idx) 
    # target_id, z = load_info(idx)
    import time
    t1 = time.time()
    print('Run estimate')
    esti_smass(ID = str(idx), mags_dict = mag_result, z = z, flag = 1, 
                if_run_gsf=True, band_as_upper = band_as_upper, metallicity=0.0,
                mag_err=[0.2]*len(mag_result), just_run = False)
    t2 = time.time()
    
    
    print(idx)
    steller_file = glob.glob('esti_smass/'+str(idx)+'/SFH_*.fits')[0]
    hdul = pyfits.open(steller_file)
    info = hdul[0].header 
    
    
    

# # print(name_list[idx])
# print( "[{0:.1f}$-${1:.1f}]".format(float(info['Mstel_16']), float(info['Mstel_84'])) )
# print( "[{0:.1f}$-${1:.1f}]".format(10**float(info['SFR_16']), 10**float(info['SFR_84'])) )
# print( "[{0:.1f}$-${1:.1f}]".format(10**float(info['T_MW_16']), 10**float(info['T_MW_84'])) )
    
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