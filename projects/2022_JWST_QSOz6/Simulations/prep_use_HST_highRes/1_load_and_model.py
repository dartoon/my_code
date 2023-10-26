#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 16:45:23 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

filt = 'f150w'
folder = 'JWST_CEERS/'
file = 'ceers5_{filt}_i2d.fits'.format(filt=filt)

im = pyfits.open(folder+file)
data_sb = im[1].data
header = im[1].header

data_sb = data_sb[:2000,:2000]
data = data_sb/header['PHOTMJSR']
value_unit = header['BUNIT']
print('For flux value in unit of MJy/sr, the zp is 31.4.') #https://en.wikipedia.org/wiki/Jansky#AB_magnitude

print("Conversion factor from {units} to DN/S for filter {f}:".format(units=header['BUNIT'], f=filt), 
      header['PHOTMJSR'])

header0 = im[0].header
img_filter = header0['FILTER']
img_cam = header0['APERNAME'] #In JDAT'simulation it is 'DETECTOR'
exptime = header0['TEXPTIME'] #The assumed exp time.
from galight.tools.astro_tools import plt_fits
plt_fits(data)

# #%% Model using galight
# from galight.data_process import DataProcess
# zp = 31.4 - 2.5*np.log10(header['PHOTMJSR'])  #Calculate the correspondingly zp as DN/S #This is wrong! See 3_...
# data_process = DataProcess(fov_image = data, target_pos = [1170., 940.], pos_type = 'pixel', header = header,
#                           rm_bkglight = False, exptime = np.ones_like(data)*exptime, if_plot=False, zp = zp)  #Gain value assuming as 1
# data_process.generate_target_materials(radius=65, create_mask = False, nsigma=2.8, if_select_obj=False,
#                                       exp_sz= 1.2, npixels = 15, if_plot=True)
# data_process.find_PSF(radius = 30, user_option = True)
# data_process.plot_overview(label = 'Example', target_label = None)

# #Start to produce the class and params for lens fitting.
# from galight.fitting_specify import FittingSpecify

# # data_process.apertures = []
# fit_sepc = FittingSpecify(data_process)
# fit_sepc.prepare_fitting_seq(point_source_num = 1) #, fix_n_list= [[0,4],[1,1]])
# fit_sepc.build_fitting_seq()

# #Plot the initial settings for fittings. 
# fit_sepc.plot_fitting_sets()

# #Setting the fitting method and run.
# from galight.fitting_process import FittingProcess
# fit_run = FittingProcess(fit_sepc, savename = 'savename', fitting_level='norm')
# fit_run.run(algorithm_list = ['PSO', 'PSO'])
# fit_run.plot_final_galaxy_fit()
# print(fit_run.final_result_galaxy[0])
# print(fit_run.final_result_ps[0])