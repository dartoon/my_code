#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 16:06:05 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import warnings, pickle
warnings.filterwarnings("ignore")
from astropy.wcs import WCS
from galight.tools.astro_tools import plt_many_fits
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import pickle, copy, glob

#%%

filt_i = 2
filt = ['F115W', 'F150W','F277W', 'F444W'][filt_i]
filefolder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/' #Folder of the target
filename = 'mosaic_nircam_f{0}w_COSMOS-Web_30mas_v0_1_i2d.fits'.format(filt[1:-1])
fitsFile = pyfits.open(filefolder+filename)
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
img = fitsFile[1].data #
wcs = WCS(header)

#%%
# RA, Dec, target_id = 149.877874999999989,  2.297638888888889, 'FMOS_J095930.7p021752'
RA, Dec, target_id = 149.933708333333300,  2.331555555555555, 'FMOS_J095944.1p021954'
# RA, Dec, target_id = 149.889708333333289,  2.358777777777778, 'FMOS_J095933.5p022132'
# RA, Dec, target_id = 149.928999999999974,  2.362527777777778, 'FMOS_J095943.0p022145'
# RA, Dec, target_id = 149.999249999999961,  2.452083333333333, 'FMOS_J095959.8p022708'

#%%
flux_mjsr = header['PHOTMJSR']
pixscale = read_pixel_scale(header)
zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
fov_noise_map = fitsFile[2].data
PSF_name = glob.glob('PSFs_library/PSF_{0}*fits'.format(filt))[0]
PSF = pyfits.getdata(PSF_name)

#%%    
from galight.tools.astro_tools import plt_fits
data_process = DataProcess(fov_image = img, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
                          rm_bkglight = False, if_plot=False, zp = zp, fov_noise_map = fov_noise_map)
data_process.generate_target_materials(radius=75, 
                                       create_mask = False, nsigma=3, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 100, if_plot=True)
data_process.noise_map[data_process.noise_map < 0.00025] = np.max(data_process.noise_map[data_process.noise_map!=0])
# import copy
# aper = data_process.apertures[0]
# aper_0 = copy.deepcopy(aper)
# aper_0.a, aper_0.b =aper_0.a/2, aper_0.b/2
# aper_1 = copy.deepcopy(aper)
# aper_1.a, aper_1.b =aper_1.a*1.5, aper_0.b*1.5
# data_process.apertures = [aper_0, aper_1]

data_process.filt = filt
del data_process.fov_image
del data_process.fov_noise_map
data_process.PSF_list = [PSF]
fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 0, supersampling_factor = 3, apertures_center_focus=True)

fit_sepc.kwargs_params['lens_light_model'][4][0]['R_sersic'] = 0.4
fit_sepc.plot_fitting_sets()
#%%

fit_run = FittingProcess(fit_sepc, savename = target_id)
fit_run.run(algorithm_list = ['PSO', 'PSO'], fitting_level=['norm','norm'])
fit_run.plot_final_galaxy_fit(target_ID =target_id, show_plot=True)
# pickle.dump(fit_run , open('fit_{1}_{0}.pkl'.format(target_id,filt), 'wb'))

# #%%Plot disk and bulge
# from galight.tools.plot_tools import total_compare
# data = fit_run.fitting_specify_class.kwargs_data['image_data']
# noise = fit_run.fitting_specify_class.kwargs_data['noise_map']
# galaxy_list = fit_run.image_host_list
# galaxy_total_image = np.zeros_like(galaxy_list[0])
# for i in range(len(galaxy_list)):
#     galaxy_total_image = galaxy_total_image+galaxy_list[i]
# model = galaxy_total_image
# norm_residual = (data - model)/noise
# flux_list_2d = [data, model, norm_residual]
# label_list_2d = ['data', 'model', 'normalized residual']
# flux_list_1d = [data, galaxy_list[0], galaxy_list[1]]
# label_list_1d = ['data', 'bulge', 'disk']
# total_compare(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d, deltaPix = fit_run.fitting_specify_class.deltaPix,
#                       zp=fit_run.zp, if_annuli=False, arrows= False, show_plot = False,
#                       target_ID = 'target_ID')

# %%For future load:
fit_run = pickle.load(open('fit_{1}_{0}.pkl'.format(target_id,filt),'rb'))
    