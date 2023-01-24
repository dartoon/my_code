#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 16:48:36 2022

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
filt = 'F356W'
idx = 0
run_folder = '../model_z6_data_id{0}/stage3_all/'.format(idx)

#%% Load data
folder = '../NIRCam_data/*/bkg_removed'.format(idx)
# info = target_info[str(idx)]
# files = glob.glob(folder+'/*{1}*{2}*{0}*.fits'.format(filt, info['target_id'][:5], info['target_id'][-4:]))
# if len(files)==1:
#     file = files[0]
# else:
#     raise ValueError('More than one files found')
# im = pyfits.open(file)
# data = im[1].data
# header = im[1].header

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
prop= ['magnitude'][0]  #Offset is in unit of arcsec.
if prop != 'offset':
    true_value = fit_run_.final_result_galaxy[0][prop]
    true_value = -2.5*np.log10(15)+fit_run_.zp  #15 is the flux I assumed.
    
else:
    fit_run_.cal_astrometry()
    true_value = np.sqrt(np.sum(np.array(fit_run_.final_result_galaxy[0]['position_xy']) - 
                   np.array(fit_run_.final_result_ps[0]['position_xy'] ))**2) * fit_run_.fitting_specify_class.deltaPix
all_sim = glob.glob('sim_result_less_exp/sim_idx{0}_{1}_seed*CombPSF.pkl'.format(idx, filt))#+\
all_sim.sort()
seedmax =  len(all_sim)
obtain_value = [] 
for seed in range(seedmax):
    fit_files = glob.glob('sim_result_less_exp/sim_idx{1}_{2}_seed{0}B*.pkl'.format(seed,idx,filt))#+\
    # fit_files = glob.glob('sim_result_less_exp/sim_idx{1}_{2}_seed{0}Based*BasedP*.pkl'.format(seed,idx,filt))#+\
    fit_files.sort()
    fit_run_list = []
    for i in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort() 
    best_sim_run = fit_run_list[sort_Chisq[0]]
    if prop != 'offset':    
        value = best_sim_run.final_result_galaxy[0][prop]
        if value - true_value < -0.6:
            print(seed, fit_files[0])
            fit_run_check=best_sim_run
            best_sim_run = fit_run_list[sort_Chisq[1]]
        value = best_sim_run.final_result_galaxy[0][prop]
    else:
        best_sim_run.cal_astrometry()
        value = np.sqrt(np.sum(np.array(best_sim_run.final_result_galaxy[0]['position_xy']) - 
                       np.array(best_sim_run.final_result_ps[0]['position_xy'] ))**2) * best_sim_run.fitting_specify_class.deltaPix
    obtain_value.append( value )
        
obtain_value = np.array(obtain_value)
plt.figure(figsize=(11, 7))
plt.hist(np.array(obtain_value[obtain_value<26.5]))
plt.axvline(x=true_value, ls='--', linewidth=1.6, c='red', zorder = 1, label='true value')
plt.ylim([0,60])
plt.xlabel(prop+' bias (obtained $-$ true)', fontsize=27)
plt.tick_params(labelsize=20)
plt.title('Simulation Result using {0} realizations'.format(seedmax), fontsize = 29)
plt.show()
print('{0:.2f}$\pm${1:.2f}'.format(np.mean(np.array(obtain_value)-true_value), np.std(np.array(obtain_value)-true_value)) )
