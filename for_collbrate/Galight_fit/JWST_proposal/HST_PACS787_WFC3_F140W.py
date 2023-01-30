#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 16:01:17 2020

@author: Xuheng Ding

You can skip this step if the QSO stamp, noise level and the PSF is ready.
"""
#photutils in version > = 0.7.2
#astropy in version astropy-4.0.1

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits


#%%
import galight.tools.astro_tools as astro_tools
import galight.tools.astro_tools as read_pixel_scale
# wht = fitsFile[2].data # The WHT map
# pixel_s = 0.03

#%%
fitsFile = pyfits.open('PACS787_hst_wfc3_f140w_030mas_drz.fits')
fov_image = fitsFile[0].data # check the back grounp
header = fitsFile[0].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
exp =  astro_tools.read_fits_exp(fitsFile[0].header)  #Read the exposure time 
# mean_wht = exp * (0.0642/0.135)**2
# exp_map = exp * wht/mean_wht

#%%Start to use galight
from galight.data_process import DataProcess
data_process = DataProcess(fov_image = fov_image, target_pos = [617, 600], pos_type = 'pixel', header = header,
                          rm_bkglight = False, exptime = exp*np.ones_like(fov_image), if_plot=False, zp = 26.450)

data_process.generate_target_materials(radius=80, create_mask = False, nsigma=2.8, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 15, if_plot=True)
#%%
data_process.deltaPix = 0.03
import copy
aper = copy.deepcopy(data_process.apertures[6] )
aper.positions = np.array([80,60])
apertures = copy.deepcopy(data_process.apertures)
apertures = apertures[:1] + [aper] + apertures[1:]
apertures[0].a = apertures[0].a*0.8
data_process.apertures = apertures

#%%PSF works.
data_process.find_PSF(radius = 40, PSF_pos_list=[[270,935]])
# data_process.find_PSF(radius = 50, PSF_pos_list = [[ 350., 1055.], [2078., 1910.]], user_option = False)
data_process.plot_overview(label = 'Example', target_label = None)

psf = data_process.PSF_list[0]
psf[psf<0] = 0.

#%%Start to produce the class and params for lens fitting.
from galight.fitting_specify import FittingSpecify
# data_process.apertures = []

fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 0) #, fix_n_list= [[0,4],[1,1]])
# psf_error_map = np.ones_like(data_process.PSF_list[data_process.psf_id_for_fitting]) *0.01 # It is in the variance unit (std^2).
# fit_sepc.prepare_fitting_seq(point_source_num = 1, psf_error_map = psf_error_map)
fit_sepc.build_fitting_seq()

#Plot the initial settings for fittings. 
fit_sepc.plot_fitting_sets()

#%%Setting the fitting method and run.
from galight.fitting_process import FittingProcess
fit_run = FittingProcess(fit_sepc, savename = 'PACS787_F140W', fitting_level=['norm','deep'])
fit_run.run(algorithm_list = ['PSO','PSO'])
fit_run.cal_astrometry()
fit_run.final_result_galaxy[:5]
fit_run.plot_final_galaxy_fit()
fit_run.dump_result()

#%%
# Test load pkl
import pickle
picklename = 'PACS787_F140W.pkl'
fitting_run_class = pickle.load(open(picklename,'rb'))
fitting_run_class.final_result_galaxy[:5]
