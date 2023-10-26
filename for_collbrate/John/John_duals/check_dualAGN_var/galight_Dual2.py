#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 15:27:48 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import glob
import pickle
#%%Cutout using online cutout tool:
import os

rerun='s21a'
# rerun='s20a'
# rerun='s19a'
print(rerun)
ra,dec = 35.2735048,  -4.6837584
point_source_num = 0  #Number of AGN in the target.
object_id = 'fit_result/'+'noPS_twosersic'
#%%Some settings for the fitting
fitting_level='deep' #shallow, deep
bands = 'I'  #TODO: Test if any band is missing.

#%%use galight to analyze:
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
print('Fitting using GaLight... ... ...')
cut_radius = 50
for band in bands:
    fitsFile = pyfits.open(glob.glob('hscdata_web_download'+'/cutout-HSC-{0}*{1}*.fits'.format(band, rerun))[0])
    file_header0 = fitsFile[0].header
    try:
        FLUXMAG0 = file_header0['FLUXMAG0']
        zp =  2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
    except:
        zp = 27.0
    PSF_file = glob.glob('hscdata_web_download'+'/psf*{1}*-{0}*.fits'.format(band, rerun))
    if PSF_file != []:
        PSF_file = PSF_file[0]
    else:
        PSF_file = glob.glob('hscdata_web_download'+'/psf*{1}*-{0}*.fits'.format(band, 's20a'))[0]
    PSF = pyfits.getdata(PSF_file)
    data_process = DataProcess(fov_image = fitsFile[1].data, fov_noise_map = fitsFile[3].data ** 0.5, target_pos = [ra, dec],
                                pos_type = 'wcs', header = fitsFile[1].header,
                                rm_bkglight = False, if_plot=False, zp = zp)
    data_process.PSF_list = [PSF]
    data_process.generate_target_materials(radius=cut_radius, create_mask = False, nsigma=2.8,
                                          exp_sz= 1.2, npixels = 25, if_plot=False)
    import copy
    apr0 = copy.deepcopy(data_process.apertures[0])
    apr1 = copy.deepcopy(data_process.apertures[0])
    apr0.positions = np.array([48., 49.])
    apr1.positions = np.array([56., 48.])
    data_process.apertures = [apr0, apr1]
    # if rerun != 's21a':
    #     picklename = 'fit_result/'+'dual_result-band-{0}-s21a.pkl'.format(band)
    #     s21_res = pickle.load(open(picklename,'rb'))
    #     data_process.target_stamp -= s21_res.image_host_list[0]
    fit_sepc = FittingSpecify(data_process)
    fit_sepc.prepare_fitting_seq(point_source_num = point_source_num, supersampling_factor=3)
    fit_sepc.plot_fitting_sets(object_id+'_fitconfig-band-{0}-{1}.png'.format(band,rerun))
    fit_sepc.build_fitting_seq()
    # if rerun != 's19a':
    #     fit_sepc.kwargs_params['point_source_model'] = s21_res.fitting_specify_class.kwargs_params['point_source_model']
    fit_run = FittingProcess(fit_sepc, savename = object_id+'_result-band-{0}-{1}'.format(band,rerun), fitting_level=fitting_level)
    fit_run.run(algorithm_list = ['PSO','PSO'], setting_list=[None,None])
    fit_run.plot_final_qso_fit(save_plot=True, target_ID= object_id +'-'+ band )
    fit_run.dump_result()
