#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 08:46:23 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
import warnings
warnings.filterwarnings("ignore")
from photutils.aperture import aperture_photometry
from photutils.aperture import RectangularAperture
from astropy.coordinates import Angle

import sys
sys.path.insert(0, '../model_z6_data_id0/')

idx = 6
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# files = glob.glob('./*fit_material*sm*/data_process_idx{0}_*_psf*.pkl'.format(idx))
# files.sort()
# run_folder = 'stage3_all/' #!!!
run_folder = '../model_z6_data_id{0}/stage3_all/'.format(idx) #!!!
# run_folder = 'stage3_all_large*/' #!!!
z_str = str(z)

import astropy.units as u
from astropy.cosmology import LambdaCDM, FlatLambdaCDM
cosmo1 = LambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, 0.7)
arc_per_kpc = cosmo1.arcsec_per_kpc_proper(z).value
filters = ['F356W', 'F150W']
# filters = ['F356W']

mags = []

import copy, matplotlib
for top_psf_id in [0]:
    for count in range(len(filters)):
        fit_run_list = []
        # idx = idx_info
        filt = filters[count]

        PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
        # idx, filt= item
        # fit_files = glob.glob(run_folder+'*fit_material*/fit_run_withcentralMask_idx{0}_{1}_FOV*.pkl'.format(idx, filt))#+\
        if idx != 1:
            fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        else:
            fit_files = glob.glob(run_folder+'*fit_material/fit_run*_fixn1_*idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        # print('idx', idx, target_id, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]
        deltaPix = fit_run.fitting_specify_class.deltaPix
        fit_run.savename = 'figures/' + fit_run.savename+'_'+filt
        count_n = 5
        
        fit_files[sort_Chisq[top_psf_id]] 
        
        Chisq_best = chisqs[sort_Chisq[top_psf_id]]
        Chisq_last= chisqs[sort_Chisq[count_n-1]]
        inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
        weight = np.zeros(len(chisqs))
        for i in sort_Chisq[:count_n]:
            weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
        
        print('idx', idx, target_id)
        print(filt)
        print(z)
        
        prop_name = 'magnitude'
        all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('host', prop_name, ": {0:.2f}+-{1:.2f}".format(weighted_value, rms_value))
        mags.append(weighted_value)

print('\nz=', z)
mag_err = [0.25]*len(mags)
fnu = [10 ** ((mags[i]-25)/(-2.5)) for i in range(len(mags))]
fnu_up = [10 ** ((mags[i]-mag_err[i]-25)/(-2.5)) for i in range(len(mags))]
fnu_dw = [10 ** ((mags[i]+mag_err[i]-25)/(-2.5)) for i in range(len(mags))]
fnu_err = [(fnu_up[i]-fnu_dw[i])/2 for i in range(len(mags))]
print('{0} {1} {2} {3}'.format(fnu[0], fnu_err[0], fnu[1], fnu_err[1]))
    # print(fnu, fnu_err)
