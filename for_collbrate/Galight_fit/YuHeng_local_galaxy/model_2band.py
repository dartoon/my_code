#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 14:06:30 2021

@author: Dartoon
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

#%%
run = 0
if run == 0:
    fitsFile = pyfits.open('./hst_mos_0074259_acs_wfc_f850lp_long_drz.fits')
    bkg_std = 0.2831354
    expzs = 1.2
    rebin_name = 'f850lp_rebin20'
elif run == 1:
    fitsFile = pyfits.open('./hst_mos_0074259_acs_wfc_f475w_drz.fits')
    bkg_std = 0.27
    expzs = 2.5
    rebin_name = 'f475w_rebin20'
fov_image = fitsFile[1].data # check the back grounp
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']

#%%
from galight.tools.astro_tools import plt_fits
numPix = np.min(fov_image.shape)
fov_image = fov_image[:numPix, :numPix]
# plt_fits(fov_image)

def rebin(image, shape,factor=None):
    sh = shape[0],image.shape[0]//shape[0],shape[1],image.shape[1]//shape[1]
    return image.reshape(sh).mean(-1).mean(1)*factor**2
factor=20
# fov_image_exp_grid=expend_grid(fov_image)
fov_image_rebin = rebin(fov_image,(int(numPix/factor),int(numPix/factor)),factor=factor)
#%%
from galight.tools.cutout_tools import cutout
PSF = cutout(fov_image, [7377, 2199], 30)
PSF = PSF[:-1,:-1]
PSF_rebin = rebin(PSF,(int(len(PSF)/factor),int(len(PSF)/factor)),factor=factor)
# plt_fits(PSF_rebin)

pyfits.PrimaryHDU(fov_image_rebin).writeto(rebin_name+'.fits',overwrite=True)

# #%%
import galight.tools.astro_tools as astro_tools
exp =  astro_tools.read_fits_exp(fitsFile[0].header)  #Read the exposure time 
exp_map = exp * np.ones_like(fov_image_rebin) #* factor**2
from galight.data_process import DataProcess
data_process = DataProcess(fov_image = fov_image_rebin, target_pos = [318.0, 239.6], pos_type = 'pixel',
                          rm_bkglight = True, exptime = exp_map, if_plot=False, zp = 27.0)
data_process.generate_target_materials(radius=len(fov_image_rebin)/2, create_mask = False, nsigma=2.8, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 300, if_plot=False, bkg_std = bkg_std)
data_process.PSF_list = [PSF_rebin]
data_process.deltaPix = 1
data_process.noise_map[np.isnan(data_process.noise_map)] = bkg_std
from galight.tools.measure_tools import mask_obj
import copy
apertures = copy.deepcopy(data_process.apertures)
apertures[0].a, apertures[0].b = apertures[0].a*expzs,  apertures[0].b*expzs
apertures[1].a, apertures[1].b = apertures[1].a*expzs,  apertures[1].b*expzs
masks = mask_obj(data_process.target_stamp, apertures, sum_mask=True)
masks = 1-masks
data_process.target_mask = masks
# plt.imshow(masks, origin='lower')
# plt.show()

#%%
savename = rebin_name+'removebkg' #+ '.pkl'
#Setting the fitting method and run.
import glob, pickle
if glob.glob(savename+'*pkl') == []:
    # apertures = copy.deepcopy(data_process.apertures)
    # apertures = [copy.deepcopy(apertures[0])] + apertures
    # apertures[0].a, apertures[0].b = apertures[0].a*0.2,  apertures[0].b*0.2
    # data_process.apertures = apertures
    from galight.fitting_specify import FittingSpecify
    fit_sepc = FittingSpecify(data_process)
    fit_sepc.prepare_fitting_seq(point_source_num = 0) #, fix_n_list= [[0,4],[1,1]])
    fit_sepc.build_fitting_seq()
    fit_sepc.plot_fitting_sets()  #The apertures shows how the images will be modelled.
    from galight.fitting_process import FittingProcess
    fit_run = FittingProcess(fit_sepc, savename = savename, fitting_level='deep')
    fit_run.run(algorithm_list = ['PSO'], setting_list=[None])
    fit_run.dump_result()
else:
    fit_run = pickle.load(open(glob.glob(savename+'*.pkl')[0],'rb'))
    fit_run.fitting_specify_class.plot_fitting_sets()

file_header = copy.deepcopy(fitsFile[1].header)
file_header['CRPIX1'] = file_header['CRPIX1']-data_process.target_pos[0]+len(data_process.target_stamp)/2
file_header['CRPIX2'] = file_header['CRPIX2']-data_process.target_pos[1]+len(data_process.target_stamp)/2

fit_run.plot_all()
# fit_run.plot_final_galaxy_fit()
res = (fit_run.flux_2d_out['data'] - fit_run.flux_2d_out['model']) #* masks
plt_fits(res)
pyfits.PrimaryHDU(res).writeto(rebin_name+'_data-model'+'.fits',overwrite=True)

