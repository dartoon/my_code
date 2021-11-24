#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 21:45:25 2020

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

#%%
import sys
# image_id = sys.argv[1] 
# band = sys.argv[2] 
# image_RA = float(sys.argv[3]) 
# image_DEC = float(sys.argv[4]) 

image_id = '000017.88+002612.6'
band = 'I'
image_RA = 0.07452999800443649
image_DEC = 0.4368380010128021

fitsFile = pyfits.open('../example_data/HSC/QSO/{0}_HSC-{1}.fits'.format(image_id, band))

fov_image= fitsFile[1].data
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
err_data= fitsFile[3].data ** 0.5

file_header0 = fitsFile[0].header
FLUXMAG0 = file_header0['FLUXMAG0']
zp =  2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
PSF = pyfits.getdata('../example_data/HSC/QSO/{0}_HSC-{1}_psf.fits'.format(image_id, band))

#%%Start to use decomprofile
from decomprofile.data_process import DataProcess

data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [image_RA, image_DEC],
                           pos_type = 'wcs', header = header,
                          rm_bkglight = True, if_plot=False, zp = zp)

data_process.noise_map = err_data

data_process.generate_target_materials(radius=None, create_mask = False, nsigma=2.8,  #radius=None, the target size would be automated set.
                                      exp_sz= 1.2, npixels = 15, if_plot=False)

data_process.PSF_list = [PSF]

data_process.checkout() #Check if all the materials is known.

#%%Start to produce the class and params for lens fitting.
from decomprofile.fitting_specify import FittingSpeficy
fit_sepc = FittingSpeficy(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 0)
fit_sepc.build_fitting_seq()

#%%Setting the fitting method and run.
from decomprofile.fitting_process import FittingProcess
savename = image_id+'-'+band
fit_run = FittingProcess(fit_sepc, savename = savename)
fit_run.run(algorithm_list = ['PSO'], setting_list = [None])  #Only PSO, not MCMC
# fit_run.plot_all()
fit_run.plot_final_galaxy_fit(target_ID= savename, show_plot = False) 
fit_run.translate_result()
# fit_run.dump_result()  #To save result as pickle file
# print(fit_run.final_result_galaxy[0])

#%%
filename_ascii = image_id + '_result.txt'
_ascii =  open(filename_ascii,'w')
_ascii.write(str(fit_run.final_result_galaxy[0]))
_ascii.close()


# #%%
# # Test load pkl
# import pickle
# picklename = savename
# fitting_run_class = pickle.load(open(picklename,'rb'))
# print(fit_run.final_result_galaxy[0])


