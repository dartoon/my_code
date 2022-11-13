#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 15:29:43 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import sys
sys.path.insert(0,'../model_z6_data_id0/')

filters =  ['F150','F356']

folder = '/Users/Dartoon/Downloads/z6JWSTNIRcam/NIRCam_J2236_stage3_all'

filt = filters[1] #!!!

idx = 1
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# #%% Remove bkg
# filter_files= glob.glob(folder+'/*{0}*{1}*i2d.fits'.format(target_id[:5], filt))  #For NIRCam
# filter_files.sort()
# from galight.data_process import DataProcess
# from galight.tools.astro_tools import read_pixel_scale
# # filename = filter_files[1]
# for filename in filter_files:
#     print("Select for", filename.split('/')[-1])
#     # Grab the JWST provided ERR map:
#     fitsFile = pyfits.open(filename)
#     fov_image = fitsFile[1].data # check the back grounp
#     header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        
#     flux_mjsr = header['PHOTMJSR']
#     pixscale = read_pixel_scale(header)
#     zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
    
#     # fov_noise_map = fitsFile[2].data
    
#     wht = fitsFile[4].data # The WHT map
#     exp = fitsFile[0].header['EFFEXPTM']
#     gain_value = 2
#     exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
        
#     from galight.data_process import DataProcess
#     data_process = DataProcess(fov_image = fov_image, target_pos = [600,1150], pos_type = 'pixel', header = header,
#                               rm_bkglight = True, exptime = exp_map, if_plot=True, zp = zp)#, fov_noise_map = fov_noise_map)
#     fitsFile[1].data = data_process.fov_image
#     fitsFile.writeto(folder+'/bkg_removed/'+filename.split('/')[-1][:-5] + '_rmbkg'+'.fits')


#%% model
filter_files= glob.glob(folder+'/bkg_removed/*{0}*{1}*i2d_rmbkg.fits'.format(target_id[:5], filt))  #For NIRCam
filter_files.sort()
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
# filename = filter_files[1]
filename = filter_files[0]
print("Select for", filename.split('/')[-1])
# Grab the JWST provided ERR map:
fitsFile = pyfits.open(filename)
fov_image = fitsFile[1].data # check the back grounp
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    
flux_mjsr = header['PHOTMJSR']
pixscale = read_pixel_scale(header)
zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux

# fov_noise_map = fitsFile[2].data

wht = fitsFile[4].data # The WHT map
exp = fitsFile[0].header['EFFEXPTM']
gain_value = 2
exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
    
from galight.data_process import DataProcess
data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
                          rm_bkglight = False, exptime = exp_map, if_plot=False, zp = zp)#, fov_noise_map = fov_noise_map)

if '150' in filt:
    data_process.generate_target_materials(radius=20, create_mask = False, nsigma=2.8, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 15, if_plot=True, cut_kernel = 'center_gaussian')
else:
    data_process.generate_target_materials(radius=45, create_mask = False, nsigma=2.8, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 15, if_plot=True, cut_kernel = 'center_gaussian')
#%%PSF works.
data_process.find_PSF(radius = 30, user_option = True, if_filter=True, nearyby_obj_filter=True)
data_process.plot_overview()
# from galight.tools.measure_tools import stack_PSF
# epsf = stack_PSF(data_process.fov_image, data_process.PSF_pos_list, psf_size=len(data_process.PSF_list[0]))
# data_process.stack_PSF()
#%%
# data_process.apertures = []
data_process.checkout() #Check if all the materials is known.
psf = data_process.PSF_list[0]
psf[psf<0] = 0.
data_process.PSF_list[0] = psf

#Start to produce the class and params for lens fitting.
from galight.fitting_specify import FittingSpecify

fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1) #, fix_n_list= [[0,4],[1,1]])
fit_sepc.build_fitting_seq()
# fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
fit_sepc.plot_fitting_sets()

#Setting the fitting method and run.
from galight.fitting_process import FittingProcess
fit_run = FittingProcess(fit_sepc, savename = 'savename')
fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
fit_run.plot_final_qso_fit(target_ID = target_id)