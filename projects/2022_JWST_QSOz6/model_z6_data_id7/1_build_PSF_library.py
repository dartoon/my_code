#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:36:24 2022

@author: Dartoon

come form 0_remove_JWST_bkg.py
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle
from galight.tools.astro_tools import plt_fits
from galight.tools.astro_tools import plt_many_fits

filters =  ['F150','F356']

data_type = 'all'
filt = filters[0] #!!!
file_NO = 0  #!!!

idx = 7

import sys
sys.path.insert(0, '../model_z6_data_id0/')
from target_info import target_info

info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

folder = '../NIRCam_data/Apr27/bkg_removed/' 
from astropy.coordinates import SkyCoord
from astropy import units as u
pos = SkyCoord('{0} {1}'.format(RA, Dec), unit=(u.hourangle, u.deg))
target_pos = np.array([pos.ra.degree, pos.dec.degree])
#%%
filter_files= glob.glob(folder+'*{0}*{1}*_rmbkg.fits'.format(target_id[:5], filt))  #For NIRCam
filter_files.sort()

from galight.tools.cutout_tools import psf_clean
filter_files = [filter_files[file_NO]]
if data_type == 'all':
    run_folder = 'stage3_{0}/'.format(data_type)
elif data_type == 'half':
    if file_NO == 0:
        run_folder = 'stage3_first_half/'
    if file_NO == 1:
        run_folder = 'stage3_second_half/'

re_select = True
clean_up = True
#%% Select PSF library:
if re_select == True:
    for filename in filter_files:
        print("Select for", filename.split('/')[-1])
        # Grab the JWST provided ERR map:
        
        fitsFile = pyfits.open(filename)
        fov_image = fitsFile[1].data # check the back grounp
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        from galight.data_process import DataProcess
        from galight.tools.astro_tools import read_pixel_scale
            
        flux_mjsr = header['PHOTMJSR']
        pixscale = read_pixel_scale(header)
        zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
        
        # fov_noise_map = fitsFile[2].data
        
        wht = fitsFile[4].data # The WHT map
        exp = fitsFile[0].header['EFFEXPTM']
        gain_value = 2
        exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
            
        data_process = DataProcess(fov_image = fov_image, target_pos = target_pos, pos_type = 'wcs', header = header,
                                  rm_bkglight = False, exptime = exp_map, if_plot=False, zp = zp)#, fov_noise_map = fov_noise_map)
        # data_process.generate_target_materials(radius=30, create_mask = False, nsigma=2.8, if_select_obj=False,
        #                                       exp_sz= 1.2, npixels = 15, if_plot=True)
        #PSF works.
        data_process.find_PSF(radius = 75, user_option = True, if_filter=False, nearyby_obj_filter=False, FWHM_sort=True)
        data_process.plot_overview()
        PSFs = data_process.PSF_list
        # PSFs = [psf_clean(PSFs[i], print_string='clean PSF '+str(i), 
        #                   ratio_to_replace=0.01, if_plot=True) for i in range(len(PSFs))]
        # data_process.PSF_list = PSFs
        # data_process.stack_PSF()
        from galight.tools.astro_tools import plt_many_fits
        plt_many_fits(PSFs)
        
        PSF_list = data_process.PSF_list
        PSF_pos_list = data_process.PSF_pos_list
        from astropy.wcs import WCS
        wcs = WCS(header)
        # sky = [wcs.pixel_to_world(PSF_pos_list[i]) for i in range(data_process.PSF_pos_list)]
        PSF_RA_DEC_list = []
        for i in range(len(PSF_pos_list)):
            RA_DEC = wcs.all_pix2world(PSF_pos_list[i][0], PSF_pos_list[i][1],0) 
            PSF_RA_DEC_list.append( [float(RA_DEC[0]), float(RA_DEC[1])] )
            
        noselect = False
        save_name = filename.split('/')[-1].split('_1_i2d')[0]
        if noselect == True:
            pickle.dump([[], [], []], open(run_folder+'material/'+save_name+'_PSF_info.pkl', 'wb'))
        else:
            pickle.dump([PSF_list, PSF_pos_list, PSF_RA_DEC_list], open(run_folder+'material/'+save_name+'_PSF_info.pkl', 'wb'))

#%% Clean up PSF:
if clean_up == True:
    PSF_lib_files = glob.glob(run_folder+'material/*'+filt+'*_PSF_info.pkl')
    PSF_lib_files.sort()
    PSF_list = []
    PSF_RA_DEC_list = []
    PSF_from_list = []  #Check if PSF are the stars
    PSF_RA_DEC_is_QSO = []
    for i in range(len(PSF_lib_files)):
        PSF_list_, PSF_pos_list_, PSF_RA_DEC_list_ = pickle.load(open(PSF_lib_files[i],'rb'))
        PSF_list = PSF_list+PSF_list_
        PSF_RA_DEC_list = PSF_RA_DEC_list + PSF_RA_DEC_list_
        PSF_from_list = PSF_from_list+ [PSF_lib_files[i]] * len(PSF_list_)
    use_PSF_RA_DEC_list = []
    # f = open("target_info.txt","r")
    # string = f.read()
    # lines = string.split('\n')   # Split in to \n
    pixscale = pixscale
    # lines = np.loadtxt('target_info.txt', dtype='str')
    pos_targets = target_pos
    idx_is_QSO = []
    for i in range(len(PSF_list)):
        psf_target_dis = np.sqrt(np.sum((np.array(PSF_RA_DEC_list[i]) - pos_targets)**2))*3600/pixscale
        if np.min(psf_target_dis) < 10:
            idx_is_QSO.append(i)
    PSF_list_clean = []
    for i, psf in enumerate(PSF_list):
        print("Clean PSF", i)
        psf = PSF_list[i]
        psf = psf_clean(psf,if_plot=False, nsigma=3, npixels=45, ratio_to_replace=0.005,
                        if_print_fluxratio=True)
        # plt_fits(psf)
        PSF_list_clean.append(psf)
    print("Before remove candidates")
    plt_many_fits(PSF_list_clean)
    manual_rm_list = [9]
    final_rm_list =  manual_rm_list + idx_is_QSO
    final_PSF_list = [PSF_list[i] for i in range(len(PSF_list)) if i not in final_rm_list]
    final_PSF_list_clean = [PSF_list_clean[i] for i in range(len(PSF_list_clean)) if i not in final_rm_list]
    final_PSF_RA_DEC_list = [PSF_RA_DEC_list[i] for i in range(len(PSF_RA_DEC_list)) if i not in final_rm_list]
    final_PSF_RA_DEC_is_QSO = [PSF_RA_DEC_list[i] for i in range(len(PSF_RA_DEC_list)) if i in idx_is_QSO]
    final_PSF_from_file_list = [PSF_from_list[i] for i in range(len(PSF_from_list)) if i not in final_rm_list]
    plt_many_fits(final_PSF_list_clean)
pickle.dump([final_PSF_list, final_PSF_list_clean, final_PSF_RA_DEC_list, final_PSF_from_file_list],
            open(run_folder+'material/'+filt+'_PSF_Library_idx{0}.pkl'.format(idx), 'wb'))
#%%
print("After remove candidates")
PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
plt_many_fits(PSF_list_clean)
