#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 16:01:17 2020

@author: Xuheng Ding

You can skip this step if the QSO stamp, noise level and the PSF is ready.
"""
#photutils in version 0.7.2
#astropy in version astropy-4.0.1


import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from matplotlib.colors import LogNorm

#%%
fitsFile = pyfits.open('../example_data/HST/QSO/1147_final_drz.fits')
fov_image = fitsFile[1].data # check the back grounp
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
wht = fitsFile[2].data # The WHT map
import decomprofile.tools_data.astro_tools as astro_tools
exp =  astro_tools.read_fits_exp(fitsFile[0].header)  #Read the exposure time 
mean_wht = exp * (0.0642/0.135)**2
exp_map = exp * wht/mean_wht

#%%Start to use decomprofile
from decomprofile.data_process import DataProcess
data_process = DataProcess(fov_image = fov_image, target_pos = [1135, 648], header = header,
                         rm_bkglight = False, exptime = exp_map, if_plot=False)
data_process.generate_target_materials(radius=60, create_mask = False, nsigma=2.8,
                                     exp_sz= 1.2, npixels = 15, if_plot=False)

#%%PSF works.
data_process.find_PSF(radius = 30, user_option = False)
# data_process.find_PSF(radius = 50, PSF_pos_list = [[ 350., 1055.], [2078., 1910.]], user_option = False)
# data_process.plot_overview(label = None)
# data_process.profiles_compare(norm_pix = 5, if_annuli=False, y_log = False,
#                   prf_name_list = (['QSO'] + ['PSF{0}'.format(i) for i in range(len(data_process.PSF_lists))]) )

#%%Start to produce the class and params for lens fitting.
from decomprofile.fitting_specify import FittingSpeficy
fit_sepc = FittingSpeficy(data_process)
fit_sepc.build_fitting_seq()

#%%Setting the fitting method and run.
from decomprofile.fitting_run import FittingRun
fit_run = FittingRun(fit_sepc)
fit_run.run()
fit_run.model_plot()
fit_run.params_corner_plot()
fit_run.flux_corner_plot()
#%%
fit_run.dump_result()


