#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 22:14:14 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
import warnings
warnings.filterwarnings("ignore")

idx = 0
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# files = glob.glob('./*fit_material*sm*/data_process_idx{0}_*_psf*.pkl'.format(idx))
# files.sort()
run_folder = 'stage3_all/' #!!!
# run_folder = 'stage3_all_large*/' #!!!
z_str = str(z)

import astropy.units as u
from astropy.cosmology import LambdaCDM, FlatLambdaCDM
cosmo1 = LambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, 0.7)
arc_per_kpc = cosmo1.arcsec_per_kpc_proper(z).value
# filters = ['F150W', 'F356W']
filters = ['F356W']
import copy, matplotlib
for top_psf_id in [0]:
    for count in range(len(filters)):
        fit_run_list = []
        # idx = idx_info
        filt = filters[count]

        PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
        # idx, filt= item
        # fit_files = glob.glob(run_folder+'*fit_material*/fit_run_withcentralMask_idx{0}_{1}_FOV*.pkl'.format(idx, filt))#+\
        fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        print('idx', idx, target_id, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]
        deltaPix = fit_run.fitting_specify_class.deltaPix
        fit_run.savename = 'figures/' + fit_run.savename+'_'+filt
        count_n = 5
        
        fit_files[sort_Chisq[top_psf_id]] 
        # fit_run.plt_
        
        Chisq_best = chisqs[sort_Chisq[top_psf_id]]
        Chisq_last= chisqs[sort_Chisq[count_n-1]]
        inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
        weight = np.zeros(len(chisqs))
        for i in sort_Chisq[:count_n]:
            weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
        
        prop_name = 'magnitude'
        # all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
        all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('host', prop_name, "{0:.2f}$\pm${1:.2f}".format(weighted_value, rms_value))
        
        prop_name = 'R_sersic'
        # all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
        all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('host Reff " ' , "{0:.2f}$\pm${1:.2f}".format(weighted_value, rms_value))
        
        print('host Reff kpc ' , "{0:.2f}$\pm${1:.2f}".format(weighted_value /arc_per_kpc, rms_value  /arc_per_kpc))
        
        prop_name = 'n_sersic'
        # all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
        all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('host n_sersic " ' , "{0:.2f}$\pm${1:.2f}".format(weighted_value, rms_value))

        prop_name = 'flux_within_frame'
        all_values = [100*(fit_run_list[i].final_result_galaxy[0][prop_name]/ (fit_run_list[i].final_result_galaxy[0][prop_name] + fit_run_list[i].final_result_ps[0][prop_name]))
                      for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('host flux ratio " ' , "{0:.1f}\%$\pm${1:.1f}\%".format(weighted_value, rms_value))

        prop_name = 'magnitude'
        all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('quasar', prop_name, "{0:.2f}$\pm${1:.2f}".format(weighted_value, rms_value))
        
        prop_name = 'flux_within_frame'
        all_values = [-2.5*np.log10(fit_run_list[i].final_result_ps[0][prop_name] + fit_run_list[i].final_result_galaxy[0][prop_name])+fit_run_list[i].zp  for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('total magnitude', prop_name, "{0:.3f}$\pm${1:.3f}".format(weighted_value, rms_value))
        
        
        prop_name = 'q'
        all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('host', prop_name, "{0:.2f}$\pm${1:.2f}".format(weighted_value, rms_value))
        
        header = fit_run.fitting_specify_class.data_process_class.header
        from astropy.wcs import WCS
        wcs = WCS(header)
        ra0, dec0 = 100, 100
        res_ra0, res_dec0 = wcs.all_pix2world([(100, 100)], 1)[0]
        xx_, yy_ = ra0, dec0
        xx_dec, yy_dec = wcs.all_world2pix([[res_ra0, res_dec0+2/3600]], 1)[0]
        # print((xx_dec - xx_), (yy_dec - yy_))
        N_angle = np.arctan((xx_dec - xx_)/(yy_dec - yy_)) * 180/np.pi
        prop_name = 'phi_G'
        all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] * 180/np.pi for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('host', prop_name, "{0:.1f}$\pm${1:.1f}".format(90-N_angle+weighted_value, rms_value))

        if 'large' in run_folder:
            for ii in [1,2]:
                prop_name = 'magnitude'
                # all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
                all_values = [fit_run_list[i].final_result_galaxy[ii][prop_name] for i in range(len(fit_run_list))]
                weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
                rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
                print('obj'+str(ii), prop_name, "{0:.2f}$\pm${1:.2f}".format(weighted_value, rms_value))
            
        
        dis = []
        for i in range(len(fit_run_list)):
            _fit_run = fit_run_list[i]
            _fit_run.cal_astrometry()
            dis.append( np.sqrt(np.sum(np.array(_fit_run.final_result_galaxy[0]['position_xy']) - 
                           np.array(_fit_run.final_result_ps[0]['position_xy'] ))**2) * deltaPix)
        all_values = np.array(dis)
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        # print('Position Offset:', round(dis[0],3), 'pixel, ', round(dis[0],2), 'kpc')
        print('positional offset (")', "{0:.2f}$\pm${1:.2f}".format(weighted_value, rms_value))

        dis = []
        for i in range(len(fit_run_list)):
            _fit_run = fit_run_list[i]
            _fit_run.cal_astrometry()
            dis.append( np.sqrt(np.sum(np.array(_fit_run.final_result_galaxy[0]['position_xy']) - 
                           np.array(_fit_run.final_result_ps[0]['position_xy'] ))**2) * deltaPix /arc_per_kpc)
        all_values = np.array(dis)
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        # print('Position Offset:', round(dis[0],3), 'pixel, ', round(dis[0],2), 'kpc')
        print('positional offset (kpc)', "{0:.2f}$\pm${1:.2f}".format(weighted_value, rms_value))

# data_process = fit_run.fitting_specify_class.data_process_class
# print(-2.5*np.log10(data_process.tbl['kron_flux'][data_process.tbl['label']==0]) + data_process.zp)