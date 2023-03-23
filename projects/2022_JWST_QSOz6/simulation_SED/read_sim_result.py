#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 17:01:12 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
folder = 'first_run/'
result_folders = glob.glob(folder + 'seed*_result')
smass_match = []
ages = []
Avs = []
for result_folder in result_folders:
    steller_file = glob.glob(result_folder.replace('result', 'sim')+'/SFH_*.fits')[0]
    hdul = pyfits.open(steller_file)
    info1 = hdul[0].header 
    smass_True = float(info1['Mstel_50'])
    ages.append(10**float(info1['T_MW_50']))
    Avs.append(float(info1['AV_50']))
    
    steller_file = glob.glob(result_folder+'/SFH_*.fits')[0]
    hdul = pyfits.open(steller_file)
    info2 = hdul[0].header 
    smass_infer = float(info2['Mstel_50'])
    # print(float(smass_infer)-float(smass_True))
    smass_match.append(float(smass_infer)-float(smass_True))
    
#%%
plt.figure(figsize=(7,5))
plt.hist(smass_match)
plt.show()

plt.figure(figsize=(7,5))
plt.scatter(ages,smass_match)
plt.show()

plt.figure(figsize=(7,5))
plt.scatter(Avs,smass_match)
plt.show()
