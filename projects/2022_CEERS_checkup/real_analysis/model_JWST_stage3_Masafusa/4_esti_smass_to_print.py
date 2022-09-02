#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 13:23:13 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from functions_for_result import esti_smass, load_prop, load_info
import glob

# ID, mags, z = 'idx0', 
# 1,2,0,51,35
# folder = '202209' #Consider only JWST
folder = '20220901' #Consider a 0.4 mag error
# folder = 'F814W_upper20220901' #Consider a 0.4 mag error
rerun = False
# idx = [1,2,0,51,35]

for idx in [35,0,2,51,1]:
# for idx in [1]:
    target_id, z = load_info(idx)
    if rerun == True:
        root_folder = '../*/*'  #Include HST
        # root_folder = './*'
        mag_result = load_prop(idx, root_folder = root_folder)
        if idx ==1 and 'F814W' in mag_result.keys():
            mag_result['F814W'] += -2.5*np.log10(5) 
            del mag_result['F606W']
            del mag_result['F125W']
            del mag_result['F160W']
        import time
        t1 = time.time()
        print('Run estimate')
        band_as_upper = []
        if idx ==35:
            band_as_upper = ['F115W']
        esti_smass(ID = folder+str(idx), mags_dict = mag_result, z = z, flag = 1, 
                    if_run_gsf=True, band_as_upper = band_as_upper, mag_err=[0.4]*len(mag_result))
        t2 = time.time()


    print(idx)
    steller_file = glob.glob('esti_smass/'+folder+str(idx)+'/SFH_*.fits')[0]
    hdul = pyfits.open(steller_file)
    info = hdul[0].header 
    print(target_id)
    print('redshift', float(info['ZMC_50']))
    print('smass', float(info['Mstel_50']) )
    print('sfr', round(10**float(info['SFR_50']),3) )
    print('age', round(10**float(info['T_MW_50']),3) )
    
    print('AV', float(info['AV_50']) )
    print('SSFR', round(float(info['SFR_50']) - float(info['Mstel_50']) ,3) )
    print("open",'esti_smass/'+folder+str(idx), '/' )
    print('\n')


# # 
# steller_file = glob.glob('esti_smass/'+folder+str(idx)+'/gsf_spec_*.fits')[0]
# hdul = pyfits.open(steller_file)
# info = hdul[0].header