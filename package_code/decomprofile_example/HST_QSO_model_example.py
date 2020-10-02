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
fitsFile = pyfits.open('../example_data/HST/QSO/1147_final_drz.fits')
fov_image = fitsFile[1].data # check the back grounp
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']

#%%
wht = fitsFile[2].data # The WHT map
import decomprofile.tools.astro_tools as astro_tools
exp =  astro_tools.read_fits_exp(fitsFile[0].header)  #Read the exposure time 
mean_wht = exp * (0.0642/0.135)**2
exp_map = exp * wht/mean_wht

#%%Start to use decomprofile
from decomprofile.data_process import DataProcess
data_process = DataProcess(fov_image = fov_image, target_pos = [1135, 648], header = header,
                          rm_bkglight = False, exptime = exp_map, if_plot=False, zp = 27.0)
data_process.generate_target_materials(radius=30, create_mask = False, nsigma=2.8,
                                      exp_sz= 1.2, npixels = 15, if_plot=False)

#%%PSF works.
data_process.find_PSF(radius = 30, user_option = False)
# data_process.find_PSF(radius = 50, PSF_pos_list = [[ 350., 1055.], [2078., 1910.]], user_option = False)
data_process.plot_overview(label = 'Example', target_label = None)
# data_process.profiles_compare(norm_pix = 5, if_annuli=False, y_log = False,
#                   prf_name_list = (['target'] + ['PSF{0}'.format(i) for i in range(len(data_process.PSF_list))]) )
data_process.checkout() #Check if all the materials is known.

#%%Start to produce the class and params for lens fitting.
from decomprofile.fitting_specify import FittingSpeficy
fit_sepc = FittingSpeficy(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1, fix_n_list= [[0,4]], fix_center_list = [[0,0]])

#Plot the initial settings for fittings. 
fit_sepc.plot_fitting_sets()

fit_sepc.build_fitting_seq()

#%%Setting the fitting method and run.
from decomprofile.fitting_process import FittingProcess
fit_run = FittingProcess(fit_sepc, savename = 'test_fix_n_4')
fit_run.run()
fit_run.plot_all()
fit_run.dump_result()
#print(fit_run.final_galaxy_result[0])
#
##%%
## Test load pkl
#import pickle
#picklename = 'result.pkl'
#fitting_run_class = pickle.load(open(picklename,'rb'))
#fitting_run_class.run_diag()

