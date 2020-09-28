#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 14:45:16 2020

@author: Xuheng Ding
"""

#photutils in version > = 0.7.2
#astropy in version astropy-4.0.1

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

#%%
fitsFile = pyfits.open('../example_data/HSC/QSO/000017.88+002612.6_HSC-I.fits')

fov_image= fitsFile[1].data
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
err_data= fitsFile[3].data ** 0.5

file_header0 = fitsFile[0].header
FLUXMAG0 = file_header0['FLUXMAG0']
zp =  2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
PSF = pyfits.getdata('../example_data/HSC/QSO/000017.88+002612.6_HSC-I_psf.fits')


#%%Start to use decomprofile
from decomprofile.data_process import DataProcess
QSO_RA = 0.07452999800443649
QSO_DEC = 0.4368380010128021
data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [QSO_RA, QSO_DEC],
                           pos_type = 'wcs', header = header,
                          rm_bkglight = False, if_plot=False, zp = zp)

data_process.noise_map = err_data

data_process.generate_target_materials(radius=35, create_mask = True, nsigma=2.8,
                                      exp_sz= 1.2, npixels = 15, if_plot=False)

data_process.PSF_list = [PSF]

data_process.checkout() #Check if all the materials is known.

#%%Start to produce the class and params for lens fitting.
from decomprofile.fitting_specify import FittingSpeficy
fit_sepc = FittingSpeficy(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
fit_sepc.build_fitting_seq()

#%%Setting the fitting method and run.
from decomprofile.fitting_process import FittingProcess
fit_run = FittingProcess(fit_sepc, savename = 'test_fix_n_4_HSC')
fit_run.run()
fit_run.plot_all()
fit_run.dump_result()
# # print(fit_run.final_galaxy_result[0])