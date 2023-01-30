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
fitsFile = pyfits.open('DEIMOS_COSMOS_761315_f277w-clear_5arcsec.fits')
fov_image = fitsFile[0].data # check the back grounp
header = fitsFile[0].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
exp =  1500 #astro_tools.read_fits_exp(fitsFile[0].header)  #Read the exposure time 
# mean_wht = exp * (0.0642/0.135)**2
# exp_map = exp * wht/mean_wht

flux_mjsr = header['PHOTMJSR']
pixscale = 0.04
zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux

#%%Start to use galight
from galight.data_process import DataProcess
data_process = DataProcess(fov_image = fov_image, target_pos = [60, 60], pos_type = 'pixel', header = header,
                          rm_bkglight = False, exptime = exp*np.ones_like(fov_image), if_plot=False, zp = zp)

data_process.generate_target_materials(radius=40, create_mask = False, nsigma=1.5, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 15, if_plot=True)
#%%
import copy
aper = copy.deepcopy(data_process.apertures[2] )
positions = [[58,41],[48,48], [40,43], [36, 38],[34,33], [27, 31],[20,30] ]
apertures = []
for i,pos in enumerate(positions):
    apre_ = copy.deepcopy(aper)
    apre_.positions = np.array(pos)
    apertures.append(apre_)
# aper.positions = np.array([80,60])
# apertures = copy.deepcopy(data_process.apertures)
# apertures = apertures[:1] + [aper] + apertures[1:]
# apertures[0].a = apertures[0].a*0.8
data_process.apertures = apertures
psf = pyfits.open('Star1_f277w-clear_2.5arcsec.fits')[0].data
psf[psf<0] = 0.
data_process.PSF_list = [psf]

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
fit_run = FittingProcess(fit_sepc, savename = 'COSMOS_761315_F277W', fitting_level=['norm','deep'])
fit_run.run(algorithm_list = ['PSO','PSO'])
fit_run.cal_astrometry()
fit_run.final_result_galaxy[:5]
fit_run.plot_final_galaxy_fit()
fit_run.dump_result()

#%%
# Test load pkl
import pickle
picklename = 'COSMOS_761315_F277W.pkl'
fitting_run_class = pickle.load(open(picklename,'rb'))
fitting_run_class.final_result_galaxy
