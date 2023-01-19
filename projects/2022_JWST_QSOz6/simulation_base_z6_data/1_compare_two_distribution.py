#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 19:14:16 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
import sys
import copy
sys.path.insert(0,'../model_z6_data_id0/')
from target_info import target_info

#%%Input number
filt = ['F356W','F150W'][1]
idx = 1
run_folder = '../model_z6_data_id{0}/stage3_all/'.format(idx)
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']


#%% Load data
folder = '../NIRCam_data/*/bkg_removed'.format(idx)

#Load PSF information:
PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))

#Load top PSF result from the overall fittings
if idx == 0:
    fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
else:
    fit_files = glob.glob(run_folder+'*fit_material*/fit_run_*idx{0}_{1}_*.pkl'.format(idx, filt))#+\
    
fit_files.sort()
fit_run_list = []
for i in range(len(fit_files)):
    fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
sort_Chisq = chisqs.argsort() 
fit_run_ = fit_run_list[sort_Chisq[0]]  # use top PSF result to run the simulation.

#%% Read Best fitting result:
prop= ['magnitude', 'R_sersic', 'n_sersic', 'offset'][0]  #Offset is in unit of arcsec.
if prop != 'offset':
    true_value = fit_run_.final_result_galaxy[0][prop]
else:
    fit_run_.cal_astrometry()
    true_value = np.sqrt(np.sum(np.array(fit_run_.final_result_galaxy[0]['position_xy']) - 
                   np.array(fit_run_.final_result_ps[0]['position_xy'] ))**2) * fit_run_.fitting_specify_class.deltaPix
all_sim = glob.glob('sim_result/sim_idx{0}_{1}_seed*CombPSF.pkl'.format(idx, filt))#+\
all_sim.sort()
seedmax =  len(all_sim)
obtain_value = [] 
obtain_value_nohost = [] 

for seed in range(seedmax):
    fit_files = glob.glob('sim_result/sim_idx{1}_{2}_seed{0}B*.pkl'.format(seed,idx,filt))#+\
    fit_files.sort()
    fit_run_list = []
    for i in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort() 
    best_sim_run = fit_run_list[sort_Chisq[0]]
    # best_sim_run.plot_final_qso_fit()
    if prop != 'offset':    
        value = best_sim_run.final_result_galaxy[0][prop]
    else:
        best_sim_run.cal_astrometry()
        value = np.sqrt(np.sum(np.array(best_sim_run.final_result_galaxy[0]['position_xy']) - 
                       np.array(best_sim_run.final_result_ps[0]['position_xy'] ))**2) * best_sim_run.fitting_specify_class.deltaPix
    obtain_value.append( value )
    
    # fit_files = glob.glob('sim_result_nohost/sim_idx{1}_{2}_seed{0}Based*.pkl'.format(seed,idx,filt))#+\
    fit_files = glob.glob('sim_result_nohost/sim_idx{1}_{2}_seed{0}Based*BasedP*.pkl'.format(seed,idx,filt))#+\
    fit_files.sort()
    fit_run_list = []
    for i in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort() 
    best_sim_run = fit_run_list[sort_Chisq[0]]
    if prop != 'offset':    
        value = best_sim_run.final_result_galaxy[0][prop]
        # print(value)
        # if value < 26:
        #     print(fit_files[sort_Chisq[0]])
        #     # if 'BasedPSF0_' not in fit_files[0]:
        #     best_sim_run.plot_final_qso_fit()
        #     input('cont')
                # fit_run_list[sort_Chisq[1]].plot_final_qso_fit()
    else:
        best_sim_run.cal_astrometry()
        value = np.sqrt(np.sum(np.array(best_sim_run.final_result_galaxy[0]['position_xy']) - 
                       np.array(best_sim_run.final_result_ps[0]['position_xy'] ))**2) * best_sim_run.fitting_specify_class.deltaPix
    if value < 25.5:
        best_sim_run = fit_run_list[sort_Chisq[1]]
        # print(seed, fit_files[0])
        fit_run = best_sim_run
        if prop != 'offset':    
            value = best_sim_run.final_result_galaxy[0][prop]
        else:
            best_sim_run.cal_astrometry()
            value = np.sqrt(np.sum(np.array(best_sim_run.final_result_galaxy[0]['position_xy']) - 
                            np.array(best_sim_run.final_result_ps[0]['position_xy'] ))**2) * best_sim_run.fitting_specify_class.deltaPix
    # if 'BasedPSF0_' not in fit_files[0]:
    obtain_value_nohost.append( value )

#%%
plt.figure(figsize=(11, 7))
plt.hist(np.array(obtain_value), color = 'orange',label='host galaxy added') #-true_value)
plt.hist(np.array(obtain_value_nohost),color='wheat', label='host not added')
plt.axvline(x=true_value, ls='--', linewidth=1.6, c='red', zorder = 1, label='true value')
# import seaborn as sns
# sns.displot(obtain_value_nohost, hist=False, color= 'red')
# sns.displot(obtain_value, hist = False, color = 'blue')


plt.xlabel('inferred '+prop, fontsize=27)
plt.tick_params(labelsize=20)
plt.title('Simulation Result for {1} ({2}) with {0} realizations'.format(seedmax, target_id, filt), fontsize = 25)
plt.legend(prop={'size':17})
plt.savefig('sim_result_two_distribution.pdf')
plt.show()
print(np.mean(np.array(obtain_value)-true_value))
print(np.std(np.array(obtain_value)-true_value))