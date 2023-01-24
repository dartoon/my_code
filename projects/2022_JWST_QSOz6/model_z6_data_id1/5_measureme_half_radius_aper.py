#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:58:31 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.insert(0,'../model_z6_data_id0/')

from galight.tools.measure_tools import flux_profile
import astropy.units as u
from astropy.cosmology import LambdaCDM
cosmo1 = LambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, 0.7)

def return_Reff_aper(fit_run, cosmo1=cosmo1):
    fit_run.cal_astrometry()
    try:
        img = fit_run.flux_2d_out['data-point source'] 
    except:
        img = fit_run.flux_2d_out['data-Point Source'] 
    if len(fit_run.image_host_list) ==1:
       img = img
    elif len(fit_run.image_host_list) ==2:
       img = img - fit_run.image_host_list[1]
    else:
        raise ValueError("Other nearby obj not subt.") 
    center =  center =  np.array([len(img)/2]*2) + np.array(fit_run.final_result_galaxy[0]['position_xy'])
    radius = len(img)/2
    r_flux, r_grids, _  =  flux_profile(img, center = center, radius = radius, 
                                        x_gridspace = None, start_p=1,
                                        q =fit_run.final_result_galaxy[0]['q'],
                                        theta = - fit_run.final_result_galaxy[0]['phi_G'],
                                        if_plot=False, fits_plot = False, grids=40)
    Reff = r_grids[r_flux>r_flux.max()/2][0] * fit_run.fitting_specify_class.deltaPix
    arc_per_kpc = cosmo1.arcsec_per_kpc_proper(z).value
    Reff_kpc = Reff/arc_per_kpc
    Reff_cir_kpc = Reff_kpc * np.sqrt(fit_run.final_result_galaxy[0]['q'])
    return Reff_kpc, Reff_cir_kpc

idx = 1
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# files = glob.glob('./*fit_material*sm*/data_process_idx{0}_*_psf*.pkl'.format(idx))
# files.sort()
run_folder = 'stage3_all/' #!!!
# run_folder = 'stage3_second_half/' #!!!
z_str = str(z)

# filters = ['F150W', 'F356W']
filters = ['F356W']
import copy, matplotlib
for top_psf_id in [0]:
    for count in range(len(filters)):
        fit_run_list = []
        # idx = idx_info
        filt = filters[count]
        
        if filt == 'F150W' :
            cmap = 'inferno'
        else:
            cmap = 'gist_heat'
        my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap)) # copy the default cmap
        my_cmap.set_bad('black')
        fit_files = glob.glob(run_folder+'fit_material/fit_run_fixn1__idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        # print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]
        # print(fit_run.final_result_galaxy)
        
        count_n = 5
        Chisq_best = chisqs[sort_Chisq[top_psf_id]]
        Chisq_last= chisqs[sort_Chisq[count_n-1]]
        inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
        weight = np.zeros(len(chisqs))
        for i in sort_Chisq[:count_n]:
            weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
        
        prop_name = 'Reff maj'
        all_values = [return_Reff_aper(fit_run_list[i])[0] for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('host', prop_name, "{0:.2f} +- {1:.2f}".format(weighted_value, rms_value))
        

        prop_name = 'Reff circ'
        all_values = [return_Reff_aper(fit_run_list[i])[1] for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('host', prop_name, "{0:.2f} +- {1:.2f}".format(weighted_value, rms_value))
        

