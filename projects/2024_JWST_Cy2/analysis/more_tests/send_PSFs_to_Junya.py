#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 20:25:37 2024

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
#Load fit_run  and Rerun

idx = 1  #!!!
import sys
sys.path.insert(0, '../../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

z_str = str(z)

# filters = ['F150W', 'F356W']
filters = ['F356W', 'F150W', 'F115W', 'F250M', 'F444W', 'F200W', 'F300M', 'F480M']
for top_psf_id in range(1):
    for count in range(len(filters)):
        fit_run_list = []
        # idx = idx_info
        filt = filters[count]
        fit_files = glob.glob('../{1}/stage3_all/*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        if filt == 'F356W' or filt == 'F150W':
            fit_files = glob.glob('../../../2022_JWST_QSOz6/model_z6_data_id{0}/stage3_all/fit_material/fit_run_fixn1__idx1_{1}*pkl'.format(idx, filt))#+\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]
        fit_run.plot_final_qso_fit()
        
        pyfits.PrimaryHDU(fit_run.fitting_specify_class.kwargs_psf['kernel_point_source']).writeto('PSF_{0}.fits'.format(filt),overwrite=True)
        
