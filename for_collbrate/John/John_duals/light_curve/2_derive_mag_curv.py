#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 11:31:18 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle, glob

name = 'PSPS+sersic_fixpos_result'
# name = 'PSPS+sersic_fixpos_nearbyPSF_result'

bands = 'GRIZY'
# band = 'G'
for band in bands:
    files = glob.glob(name+"*band{0}*pkl".format(band))
    files.sort()
    mag_list = []
    date_list = []
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
        mag_list.append([mag0, mag1])
        date = date[2:4]+date[5:7]+date[8:10]
        date_list.append(date)
        
    mag_list = np.array(mag_list)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    # for i in range(len(date_list)):
    plt.plot(np.arange(len(date_list)), mag_list[:,0], label = 'left PS')
    plt.plot(np.arange(len(date_list)), mag_list[:,1], label = 'right PS')
    if band != 'Z' and band != 'I':
        for i in range(len(date_list)):
            plt.text(i, np.max(mag_list), date_list[i], label = 'left PS',fontsize=14)
    else:
        print("Test")
        for i in np.arange(0, len(date_list),3):
            plt.text(i, np.max(mag_list), date_list[i], label = 'left PS',fontsize=14)
    plt.title('Band '+band,fontsize=20)
    
    # labels = [item.get_text() for item in ax.get_xticklabels()]
    # labels = date_list
    # ax.set_xticklabels(labels)
    plt.tick_params(labelsize=20)
    plt.legend(loc=4,prop={'size':20})
    plt.savefig('light_curve_band'+band+'.pdf')
    plt.show()