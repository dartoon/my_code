#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 16:53:07 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob 

from galight.data_process import DataProcess
from galight.tools.astro_tools import plt_fits
from galight.tools.astro_tools import plt_many_fits
from astropy.wcs import WCS
import pickle
folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_HST_data/'
# all_files= glob.glob(folder+'/egs_all_wfc3_ir_*_030mas_v1.9_drz.fits')  #For IR
all_files= glob.glob(folder+'/egs_all_acs_wfc*_030mas_v1.9_drz.fits')  #For ACS
all_files.sort()

#%%
i = 0

fitsFile = pyfits.open(all_files[i])
fov_image = fitsFile[0].data # check the back grounp
header = fitsFile[0].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
if i == 0:
    filt = header['FILTER1']
else:
    filt = header['FILTER2']
# filt = header['filter']
#%%
# plt_fits(fov_image[14000:14000+8000,19000:19000+8000])

data_process = DataProcess(fov_image = fov_image, target_pos = [100, 100], 
                            pos_type = 'wcs', header = header,
                            rm_bkglight = False, if_plot=False, zp = 27 )

data_process.find_PSF(radius = 75, user_option = True, if_filter=False, 
                      nearyby_obj_filter=False, FWHM_sort=True)
# data_process.plot_overview()
PSFs = data_process.PSF_list
plt_many_fits(PSFs)

PSF_list = data_process.PSF_list
PSF_pos_list = data_process.PSF_pos_list
wcs = WCS(header)
# sky = [wcs.pixel_to_world(PSF_pos_list[i]) for i in range(data_process.PSF_pos_list)]
PSF_RA_DEC_list = []
save_name = filt + '_all'
for i in range(len(PSF_pos_list)):
    RA_DEC = wcs.all_pix2world(PSF_pos_list[i][0], PSF_pos_list[i][1],0) 
    PSF_RA_DEC_list.append( [float(RA_DEC[0]), float(RA_DEC[1])] )
hold = input('Hold ... OK?\n')
pickle.dump([PSF_list, PSF_pos_list, PSF_RA_DEC_list], open('material/'+save_name+'_PSF_info.pkl', 'wb'))

#%% Clean up PSF:
from galight.tools.astro_tools import read_pixel_scale
from galight.tools.cutout_tools import psf_clean
clean_up = True
if clean_up == True:
    PSF_list = []
    PSF_RA_DEC_list = []
    PSF_RA_DEC_is_QSO = []
    
    PSF_list, PSF_pos_list, PSF_RA_DEC_list = pickle.load(open('material/'+save_name+'_PSF_info.pkl','rb'))
    use_PSF_RA_DEC_list = []
    f = open("../model_JWST_stage3_Masafusa/material/target_info.txt","r")
    string = f.read()
    lines = string.split('\n')   # Split in to \n
    pixscale = read_pixel_scale(header)
    lines = np.loadtxt('../model_JWST_stage3_Masafusa/material/target_info.txt', dtype='str')
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
    # plt_many_fits(final_PSF_list_clean)
hold = input('Hold final one ... OK?\n')
pickle.dump([final_PSF_list, final_PSF_list_clean, final_PSF_RA_DEC_list],
            open('material/'+filt+'_PSF_Library.pkl', 'wb'))