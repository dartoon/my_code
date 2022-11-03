#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 09:06:09 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import sys
import pickle
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.tools.measure_tools import measure_bkg
from astropy.wcs import WCS
from galight.tools.cutout_tools import common_data_class_aperture
from galight.tools.plot_tools import plot_data_apertures_point
from galight.tools.cutout_tools import cutout
import warnings
warnings.filterwarnings("ignore")
from target_info import target_info
from galight.tools.astro_tools import plt_fits, plt_many_fits
#%%
data_type = 'all' 
filt = 'F356W'
file_NO = 0

idx = 0
folder = '/Users/Dartoon/Downloads/z6JWSTNIRcam/NIRCam_J2255_stage3_{0}/bkg_removed'.format(data_type)
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
# jwst_all_filenames = glob.glob(folder+'/*{0}*{1}*.fits'.format(target_id[:5], filts[0]))
jwst_all_filenames = glob.glob(folder+'/*{0}*.fits'.format(filt))
jwst_all_filenames.sort()
file = jwst_all_filenames[file_NO]
if data_type == 'all':
    run_folder = 'stage3_{0}/'.format(data_type)
elif data_type == 'half':
    if file_NO == 0:
        run_folder = 'stage3_first_half/'
    if file_NO == 1:
        run_folder = 'stage3_second_half/'
result_folder = run_folder + 'fit_result/'

#%%Demonstrate the local envs.
cut_kernel = None #After pos correct then, do the nearest_obj_center
# filts = ['F356W', 'F150W']
for filt in [filt]:
    _fitsFile = pyfits.open(file)
    fov_image = _fitsFile[1].data # check the back grounp
    header = _fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    flux_mjsr = header['PHOTMJSR']
    pixscale = read_pixel_scale(header)
    zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
    wht = _fitsFile[4].data # The WHT map
    exp = _fitsFile[0].header['EFFEXPTM']
    print("Exp time:", exp)
    if _fitsFile[0].header['CHANNEL'] == 'LONG':
        gain_value = 2
        expsize = 1
        exppix = 1
    else:
        gain_value = 1.8
        expsize = 1
        exppix = 2
    exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
    data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], #The final cut center is dete
                                   pos_type = 'wcs', header = header,rm_bkglight = False, 
                                   if_plot=False, zp = zp, exptime= exp_map, 
                                   fov_noise_map = None)
    
    data_process.generate_target_materials(radius=180 * expsize, create_mask = False, nsigma=1.5, 
                                            cut_kernel = None, skip = True)
    # plt_fits(data_process.target_stamp)
    data_process.apertures = []
    data_process.plot_aperture()
# input('press Enter')
#%% Show PSF:
print("Show how PSF are:")

PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
data_process.find_PSF(radius = 15, PSF_pos_list = PSF_RA_DEC_list, pos_type = 'wcs')
data_process.plot_overview()
# plt_many_fits(PSF_list_clean)
plt_many_fits(data_process.PSF_list)
# input('press Enter')
#%%
# fit_files = glob.glob('stage3_*/'+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
fit_files.sort()
fit_run_list = []
for i in range(len(fit_files)):
    fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
sort_Chisq = chisqs.argsort()  
for ii, top_psf_id in enumerate([0]):
    if ii == 0:
        fit_run_list[0].fitting_specify_class.plot_fitting_sets()
    # fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_CombPsfs*.pkl'.format(idx, filt))#+\
    print(fit_files[sort_Chisq[top_psf_id]])
    fit_run = fit_run_list[sort_Chisq[top_psf_id]]
    fit_run.savename = 'figures/' + fit_run.savename+'_'+filt
    fit_run.plot_final_qso_fit(target_ID = target_id+'$-$'+filt, save_plot = False)
    fit_run.final_result_galaxy
    
    fit_run.cal_astrometry()
    dis = np.sqrt(np.sum(np.array(fit_run.final_result_galaxy[0]['position_xy']) - 
                   np.array(fit_run.final_result_ps[0]['position_xy'] ))**2)
    import astropy.units as u
    from astropy.cosmology import LambdaCDM, FlatLambdaCDM
    cosmo1 = LambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, 0.7)
    arc_per_kpc = cosmo1.arcsec_per_kpc_proper(6.3).value
    print('Position Offset:', round(dis,3), 'pixel, ', round(dis * fit_run.fitting_specify_class.deltaPix /arc_per_kpc ,2), 'kpc')