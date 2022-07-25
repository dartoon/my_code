#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 16:13:13 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle

folder = 'output/'

# ID = '142008.61+530004.0'
# ID = 'aegis_630'
# filt = 'f356w'
# filt = 'f150w'

#Should use the use_PSF_list info to rule out the QSO.
import glob
target_id = 'aegis_511'
pkl_files = glob.glob(folder+'*fit_result_PsfLib_'+target_id+'.pkl')
pkl_files.sort()

import copy, matplotlib
# cmap = ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
#             'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
#             'hot', 'afmhot', 'gist_heat', 'copper']

cmap = ['winter','summer','afmhot','spring', 'autumn', 'gist_heat','hot' ]

for i, filename in enumerate(pkl_files):
    fit_run_list, use_PSF_list = pickle.load(open(filename,'rb'))
    filt = filename.split('output/')[1].split('_')[0]
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    idx_counts = chisqs.argsort()
    if chisqs[idx_counts[0]] < chisqs[idx_counts[1]]/2.5:
        print("The AGN org image is used PSF in", target_id , filt)
        idx_counts = idx_counts[1:]
    fit_run = fit_run_list[idx_counts[0]]
    my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap[i])) # copy the default cmap
    my_cmap.set_bad('black')
    print(cmap[i])
    fit_run.plot_final_qso_fit(target_ID = target_id + ' '+filt, cmap = my_cmap)
    host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
    AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
    ratio = host_flux/(host_flux+AGN_flux)
    print(filt, round(fit_run.final_result_galaxy[0]['flux_within_frame'],2),
          "host ratio",round(ratio,2),
          round(fit_run.final_result_galaxy[0]['magnitude'],2) )