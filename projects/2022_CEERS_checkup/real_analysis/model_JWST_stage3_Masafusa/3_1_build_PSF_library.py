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

filters =  ['F200', 'F444']

folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_JWST_Masafusa/bkg_removed'

filt = filters[1]

filter_files= glob.glob(folder+'/NIRCam2_?_*{0}*.fits'.format(filt))  #For NIRCam
filter_files.sort()

from galight.tools.cutout_tools import psf_clean

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
            
        data_process = DataProcess(fov_image = fov_image, target_pos = [0,0], pos_type = 'pixel', header = header,
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
            pickle.dump([[], [], []], open('material/'+save_name+'_PSF_info.pkl', 'wb'))
        else:
            pickle.dump([PSF_list, PSF_pos_list, PSF_RA_DEC_list], open('material/'+save_name+'_PSF_info.pkl', 'wb'))

#%% Clean up PSF:
    
filt = 'F444'
if clean_up == True:
    PSF_lib_files = glob.glob('material/*'+filt+'*_PSF_info.pkl')
    if filt == 'F200' or filt == 'F444':
        remove_i = [i for i in range(len(PSF_lib_files)) if 'NIRCam2_F' in PSF_lib_files[i]][0]
        del PSF_lib_files[remove_i]
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
    f = open("material/target_info.txt","r")
    string = f.read()
    lines = string.split('\n')   # Split in to \n
    
    # pixscale = data_process.deltaPix
    if int(filt[1:]) > 200:
        pixscale = 0.031 #!!!
    else:
        pixscale = 0.0156
    lines = np.loadtxt('material/target_info.txt', dtype='str')
    pos_targets = lines[:, 1:3].astype(np.float64)
    idx_is_QSO = []
    for i in range(len(PSF_list)):
        psf_target_dis = np.sqrt(np.sum((np.array(PSF_RA_DEC_list[i]) - pos_targets)**2, axis=1))*3600/pixscale
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
    # plt_many_fits(final_PSF_list_clean)
pickle.dump([final_PSF_list, final_PSF_list_clean, final_PSF_RA_DEC_list, final_PSF_from_file_list],
            open('material/'+filt+'_PSF_Library.pkl', 'wb'))
#%%
print("After remove candidates")
PSF_lib_files = glob.glob('material/*'+filt[:-1]+'*_PSF_Library.pkl')[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
plt_many_fits(PSF_list_clean)
