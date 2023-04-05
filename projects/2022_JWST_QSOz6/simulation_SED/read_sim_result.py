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

# path = 'first_run/'  #Av 0.3-1, Metal, -1, -0.3, No Emission line.
# path = 'second_run/'  #Av 0.3-1, Metal, -1, -0.3, logU -3, -1
# path = 'third_run/'   #Av 0.3-2, Metal, -1, -0.3
# path = 'fourth_run/'   #Av 0.3-3, Metal, -2, 0, Age, 0.01,0.75
# path = 'prior1_run/'   #Av 0.3-1, Metal, -1, 0.3, Age, 0.01,0.75
# path = 'prior2_run/'   #Av 0.3-1.5, Metal, -2, 0, Age, 0.01,0.75
# path = 'prior2_run_Av1.5/'   #Av 0.3-1.5, Metal, -2, 0, Age, 0.01,0.75
path = 'prior1_run_idx0/'   #Av 0.3-1, Metal, -2, 0, Age, 0.01,0.75

folder = path+'F150W_F356W/'  #F150W_F356W
# folder = path+'F150W_F200W/'
# folder = path+'F150W_F277W/'
# folder = path+'F277W_F356W/'
# folder = path+'F356W_F444W/'
# folder = path+'F150W_F277W_F356W/'
# folder = path+'F150W_F277W_F356W_F444W/' 

result_folders = glob.glob(folder + 'seed*_result')

sim_folder = path + 'simulation/'
smass_match = []
age_match = []
ages = []
Avs = []

sim_mass_up_low = []
mod_mass_up_low = []
for result_folder in result_folders:
    seed = result_folder.split('seed')[1].split('_')[0]
    steller_file = glob.glob(sim_folder+'seed{0}_sim/SFH_*.fits'.format(seed))[0]
    hdul = pyfits.open(steller_file)
    info1 = hdul[0].header 
    smass_True = float(info1['Mstel_50'])
    sim_mass_up_low.append(float(info1['Mstel_84']) - float(info1['Mstel_50']) )
    ages.append(10**float(info1['T_MW_50']))
    Avs.append(float(info1['AV_50']))
    
    steller_file = glob.glob(result_folder+'/SFH_*.fits')[0]
    hdul = pyfits.open(steller_file)
    info2 = hdul[0].header 
    smass_infer = float(info2['Mstel_50'])
    mod_mass_up_low.append(float(info2['Mstel_84']) - float(info2['Mstel_50']) )
    # print(float(smass_infer)-float(smass_True))
    smass_match.append(float(smass_infer)-float(smass_True))
    # if len(seed) == 1:
    # # if seed == '7':
    #     print(seed)
    #     print(float(info1['Mstel_84']) - float(info1['Mstel_16']))
    #     print(float(info2['Mstel_84']) - float(info2['Mstel_16']))
    age_match.append( 10**float(info2['T_MW_50'])-10**float(info1['T_MW_50']) )
    
#%%
Avs = np.array(Avs)
ages = np.array(ages)
smass_match = np.array(smass_match)
sim_mass_up_low = np.array(sim_mass_up_low)
mod_mass_up_low = np.array(mod_mass_up_low)
comb= folder.split('/')[-2]


plt.figure(figsize=(7,5))
if 'first' in path:
    eff_bool = (sim_mass_up_low<100)
elif 'idx0' in path:
    eff_bool = (sim_mass_up_low<100)
else:
    # eff_bool = (sim_mass_up_low<0.4) *(mod_mass_up_low<0.4) #* (sim_mass_up_low>0.15) *(mod_mass_up_low>0.15)
    eff_bool = (sim_mass_up_low<np.percentile(sim_mass_up_low,70)) *(mod_mass_up_low<np.percentile(mod_mass_up_low,70)) #* (sim_mass_up_low>0.15) *(mod_mass_up_low>0.15)
plt.hist(smass_match[eff_bool])
plt.xlabel('log(M*) scatter (infer $-$ True): {0:.2f}$\pm${1:.2f}'.format(np.mean(smass_match[eff_bool]), np.std(smass_match[eff_bool])),fontsize=15)
plt.title("filter set: "+comb, fontsize=15)
plt.show()
print(folder)
print(len(smass_match[eff_bool])/len(smass_match))
print(len(smass_match))
# plt.hist(Avs[sim_mass_up_low<0.6])
# plt.hist(sim_mass_up_low[sim_mass_up_low<np.percentile(sim_mass_up_low,60)])

# plt.figure(figsize=(7,5))
# plt.hist(age_match)
# plt.xlabel('age scatter (infer $-$ True): $\pm${0:.2f} Gyr'.format(np.std(age_match)),fontsize=15)
# plt.title("filter set: "+comb, fontsize=15)
# plt.show()
# print(folder)
# print(np.std(age_match))

# plt.figure(figsize=(7,5))
# plt.scatter(ages,smass_match)
# plt.show()

# plt.figure(figsize=(7,5))
# plt.scatter(Avs,smass_match)
# plt.show()
