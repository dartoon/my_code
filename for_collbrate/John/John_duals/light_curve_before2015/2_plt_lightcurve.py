#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 11:21:24 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from _0_scp_data import data_list as obs_list
import pickle 


mag_list = []
date_list = []
for item in obs_list:
    b ,c = item[2:4]
    file = 'HSC_data/CORR-*{0}-*{1}.fits'.format(b, c)
    name = 'HSC_{0}_{1}'.format(b, c)
    file = glob.glob(file)[0]
    # 'HSC_{0}_{1}'.format(b, c)
    pkl_file = 'results/HSC_{0}_{1}.pkl'.format(b, c)
    if glob.glob(pkl_file) != []:
        # print(pkl_file)
        fit_run = pickle.load(open(pkl_file,'rb'))
        ps_result = fit_run.final_result_ps
        mag0, mag1 = ps_result[0]['magnitude'], ps_result[1]['magnitude']
        if ps_result[0]['wcs_RaDec'][0] - ps_result[1]['wcs_RaDec'][0] > 0 and ps_result[0]['wcs_RaDec'][1] - ps_result[1]['wcs_RaDec'][1]>0:
            mag_list.append([mag0, mag1])
            date = item[0]
            date = date[2:4]+date[5:7]+date[8:10]
            date_list.append(date)
        else:
            print(file, 'mis obj')
            print('Warning, miss shift PS0 PS1', ps_result[0]['ra_image'], ps_result[1]['ra_image'])

date_list_ = np.array(date_list, dtype=np.int0)
mag_list = np.array(mag_list)
date_list_0 = [140922, 141123, 150811]
mag_av_list = []
mag_std_list = []
for date in date_list_0:
    mag_av_list.append( np.mean(mag_list[date_list_ == date], axis = 0))
    mag_std_list.append( np.std(mag_list[date_list_ == date], axis = 0))
mag_std_list = np.array(mag_std_list)
#%%
files = glob.glob("../light_curve/results/*pkl")
files.sort()
mag_list_1 = []
date_list_1 = []
for file in files:
    date = file.split('result-')[1][:10]
    filename = file
    fit_run = pickle.load(open(filename,'rb'))
    # fit_run.plot_final_qso_fit()
    ps_result = fit_run.final_result_ps
    mag0, mag1 = ps_result[0]['magnitude'], ps_result[1]['magnitude']
    
    if ps_result[0]['ra_image']<0:
        print(file)
        print('Warning, miss PS0 PS1', ps_result[0]['ra_image'], ps_result[1]['ra_image'])
        fit_run.plot_final_qso_fit()
        fit_run.fitting_specify_class.plot_fitting_sets()
    # print(date, mag0, mag1)
    # print(ps_result[0]['ra_image'], ps_result[1]['ra_image'])
    mag_list_1.append([mag0, mag1])
    date = date[2:4]+date[5:7]+date[8:10]
    date_list_1.append(int(date))

#%%

mag_list = mag_av_list + mag_list_1
date_list = date_list_0+date_list_1
mag_list = np.array(mag_list)

fig, ax = plt.subplots(figsize=(12, 6))
# for i in range(len(date_list)):
plt.plot(np.arange(len(date_list)), mag_list[:,0], label = 'left PS', c='red')
plt.plot(np.arange(len(date_list)), mag_list[:,1], label = 'right PS', c='blue')
plt.errorbar(np.arange(3), mag_list[:3,0], yerr= mag_std_list[:3,0], c='red')
plt.errorbar(np.arange(3), mag_list[:3,1], yerr= mag_std_list[:3,1], c='blue')
# if band != 'Z' and band != 'I':
for i in range(len(date_list)):
    plt.text(i-0.4, 18.8 - 0.1*(i%5), date_list[i],fontsize=14)
plt.title('Band '+'i',fontsize=20)

plt.tick_params(labelsize=20)
plt.legend(loc=3,prop={'size':20})
ax.xaxis.set_ticks([])
# plt.savefig('light_curve_band'+band+'.pdf')
plt.show()