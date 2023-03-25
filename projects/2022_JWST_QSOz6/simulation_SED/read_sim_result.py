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

path = 'second_run/'
folder = path+'F150W_F356W/'  #F150W_F356W
# folder = path+'F150W_F200W/'
# folder = path+'F150W_F277W/'
# folder = path+'F277W_F356W/'
# folder = path+'F356W_F444W/'
# folder = path+'/F150W_F277W_F356W/'
# folder = path+'/F150W_F277W_F356W_F444W/'

result_folders = glob.glob(folder + 'seed*_result')

sim_folder = 'first_run/simulation/'
smass_match = []
age_match = []
ages = []
Avs = []

check_mass_up_low = []
for result_folder in result_folders:
    seed = result_folder.split('seed')[1].split('_')[0]
    steller_file = glob.glob(sim_folder+'/seed{0}_sim/SFH_*.fits'.format(seed))[0]
    hdul = pyfits.open(steller_file)
    info1 = hdul[0].header 
    smass_True = float(info1['Mstel_50'])
    check_mass_up_low = float(info1['Mstel_50']) - float(info1['Mstel_50']) 
    ages.append(10**float(info1['T_MW_50']))
    Avs.append(float(info1['AV_50']))
    
    steller_file = glob.glob(result_folder+'/SFH_*.fits')[0]
    hdul = pyfits.open(steller_file)
    info2 = hdul[0].header 
    smass_infer = float(info2['Mstel_50'])
    # print(float(smass_infer)-float(smass_True))
    smass_match.append(float(smass_infer)-float(smass_True))
    age_match.append( 10**float(info2['T_MW_50'])-10**float(info1['T_MW_50']) )
    
#%%
comb= folder.split('/')[-2]
plt.figure(figsize=(7,5))
plt.hist(smass_match)
plt.xlabel('log(M*) scatter (infer $-$ True): $\pm${0:.2f}'.format(np.std(smass_match)),fontsize=15)
plt.title("filter set: "+comb, fontsize=15)
plt.show()
print(folder)
print(np.std(smass_match))

plt.figure(figsize=(7,5))
plt.hist(age_match)
plt.xlabel('age scatter (infer $-$ True): $\pm${0:.2f} Gyr'.format(np.std(age_match)),fontsize=15)
plt.title("filter set: "+comb, fontsize=15)
plt.show()
print(folder)
print(np.std(age_match))

# plt.figure(figsize=(7,5))
# plt.scatter(ages,smass_match)
# plt.show()

# plt.figure(figsize=(7,5))
# plt.scatter(Avs,smass_match)
# plt.show()
