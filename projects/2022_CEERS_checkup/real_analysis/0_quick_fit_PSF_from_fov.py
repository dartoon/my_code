#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 15:48:39 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o003_t023_nircam_clear-f356w/jw01345-o003_t023_nircam_clear-f356w_i2d.fits
# aegis_463 214.7768 52.825876 flux: 1.0515807 redshift: 2.274 2.24
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o003_t023_nircam_clear-f356w/jw01345-o003_t023_nircam_clear-f356w_i2d.fits
# aegis_482 214.75522 52.836795 flux: 0.9378855 redshift: 3.465 3.345
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o004_t024_nircam_clear-f356w/jw01345-o004_t024_nircam_clear-f356w_i2d.fits
# aegis_477 214.87073 52.833117 flux: 0.5838029 redshift: 2.317 2.136
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-c1000_t021_nircam_clear-f356w/jw01345-c1000_t021_nircam_clear-f356w_i2d.fits
# 142008.61+530004.0 215.03590035000232 53.00111935271983 0.64892083 2.588


# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o003_t023_nircam_clear-f150w/jw01345-o003_t023_nircam_clear-f150w_i2d.fits
# aegis_463 214.7768 52.825876 flux: 1.7630914 redshift: 2.274 2.24
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o003_t023_nircam_clear-f150w/jw01345-o003_t023_nircam_clear-f150w_i2d.fits
# aegis_482 214.75522 52.836795 flux: 1.5896 redshift: 3.465 3.345
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o004_t024_nircam_clear-f150w/jw01345-o004_t024_nircam_clear-f150w_i2d.fits
# aegis_477 214.87073 52.833117 flux: 0.28244677 redshift: 2.317 2.136
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-c1001_t021_nircam_clear-f150w/jw01345-c1001_t021_nircam_clear-f150w_i2d.fits
# 142008.61+530004.0 215.03590035000232 53.00111935271983 0.23676294 2.588

#%%
pick = 0
filt = 150
if pick == 0:
    if filt != 150:
        fitsFile = pyfits.open('/Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o003_t023_nircam_clear-f356w/jw01345-o003_t023_nircam_clear-f356w_i2d.fits')
    elif filt == 150:
        fitsFile = pyfits.open('/Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o003_t023_nircam_clear-f150w/jw01345-o003_t023_nircam_clear-f150w_i2d.fits')
    target_id = 'aegis_463'
    RA, Dec = 214.7768, 52.825876
elif pick ==1:
    if filt != 150:
        fitsFile = pyfits.open('/Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o003_t023_nircam_clear-f356w/jw01345-o003_t023_nircam_clear-f356w_i2d.fits')
    elif filt == 150:
        fitsFile = pyfits.open('/Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o003_t023_nircam_clear-f150w/jw01345-o003_t023_nircam_clear-f150w_i2d.fits')
    target_id = 'aegis_482'
    RA, Dec = 214.75522, 52.836795
elif pick ==2:
    if filt != 150:
        fitsFile = pyfits.open('/Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o004_t024_nircam_clear-f356w/jw01345-o004_t024_nircam_clear-f356w_i2d.fits')
    elif filt == 150:
        fitsFile = pyfits.open('/Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o004_t024_nircam_clear-f150w/jw01345-o004_t024_nircam_clear-f150w_i2d.fits')
    target_id = 'aegis_477'
    RA, Dec = 214.87073, 52.833117 
elif pick ==3:
    if filt != 150:
        fitsFile = pyfits.open('/Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-c1000_t021_nircam_clear-f356w/jw01345-c1000_t021_nircam_clear-f356w_i2d.fits')
    elif filt == 150:
        fitsFile = pyfits.open('/Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-c1001_t021_nircam_clear-f150w/jw01345-c1001_t021_nircam_clear-f150w_i2d.fits')
    target_id = '142008.61+530004.0'
    RA, Dec = 215.03590035000232, 53.00111935271983

fov_image = fitsFile[1].data # check the back grounp
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']

#%% Grab the JWST provided ERR map:

from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
    
flux_mjsr = header['PHOTMJSR']
pixscale = read_pixel_scale(header)

zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
wht = fitsFile[4].data # The WHT map
# exp =  header['XPOSURE']  #Read the exposure time 
exp = fitsFile[0].header['EFFEXPTM']
gain_value = 2
exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
                          rm_bkglight = True, if_plot=False, zp = zp, exptime= exp_map )
data_process.generate_target_materials(radius=30, create_mask = False, nsigma=2.8, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 15, if_plot=True)
if pick ==2:
    data_process.target_mask = data_process.target_stamp != 0
    data_process.noise_map = np.nan_to_num(data_process.noise_map, nan=1000)
# 
#%%PSF works.
data_process.find_PSF(radius = 40, user_option = True, if_filter=True, neighborhood_size=50,
                      nearyby_obj_filter=True, FWHM_sort=True, select_all=False)
data_process.plot_overview(label = 'Example', target_label = None)

#%%
data_process.profiles_compare(norm_pix = 5, if_annuli=False, y_log = False, radius = 60,
                  prf_name_list = (['target'] + ['PSF{0}'.format(i) for i in range(len(data_process.PSF_list))]) )
data_process.checkout() #Check if all the materials is known.

#%%Start to produce the class and params for lens fitting.
from galight.fitting_specify import FittingSpecify

fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1) #, fix_n_list= [[0,4],[1,1]])

fit_sepc.build_fitting_seq()

from galight.fitting_process import FittingProcess
fit_run = FittingProcess(fit_sepc, savename = 'savename', fitting_level=['norm','norm'])
fit_run.run(algorithm_list = ['PSO','PSO'])
            # setting_list = [{'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}, {'n_burn': 100, 'n_run': 100, 'walkerRatio': 10,'sigma_scale': .1}])
fit_run.plot_final_qso_fit()
print(fit_run.final_result_galaxy)
print(fit_run.final_result_ps)
