#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 16:11:03 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import glob, pickle
from galight.tools.cutout_tools import psf_clean
from astropy.wcs import WCS
from galight.tools.astro_tools import plt_fits
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.tools.astro_tools import plt_many_fits

filt_i = 3
filt = ['F115W', 'F150W','F277W', 'F444W'][filt_i]
filefolder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/'
filename = 'mosaic_nircam_f{0}w_COSMOS-Web_30mas_v0_1_i2d.fits'.format(filt[1:-1])
# filename = '1727_cosmos_mosaic_miri_exptime_scale1.0.fits'
fitsFile = pyfits.open(filefolder+filename)
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
img = fitsFile[1].data #
print(img.shape)
wcs = WCS(header)

#%%
flux_mjsr = header['PHOTMJSR']
pixscale = read_pixel_scale(header)
zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
rx, ry = img.shape
cut_p = 4 #cut_p
x_size, y_size = int(rx/4+100), int(ry/4+100)
all_PSF_pos_list_ = []
FWHM_filer = [2.1,2.25,4.5, 5.5][filt_i]  #Set 2.1 as upper limit
# flux_filt = [[], []]
rerun = False
if rerun == True:
    for i in range(cut_p):
        for j in range(cut_p):
            print('x',x_size*i,x_size*(i+1), 'y', y_size*j,y_size*(j+1))
            img_i = img[x_size*i:x_size*(i+1), y_size*j:y_size*(j+1)]
            ct_pos_i, ct_pos_j =  x_size*i, y_size*j
            
            data_process = DataProcess(fov_image = img_i, target_pos = [10,10], pos_type = 'pixel', header = header,
                                      rm_bkglight = False, if_plot=False, zp = zp)#, fov_noise_map = fov_noise_map)
            # data_process.generate_target_materials(radius=30, create_mask = False, nsigma=2.8, if_select_obj=False,
            #                                       exp_sz= 1.2, npixels = 15, if_plot=True)
            #PSF works.
            if img_i.max() == 0:
                print("ignore", i,j)
                continue
            data_process.find_PSF(radius = 50, user_option = True, if_filter=True,
                                  FWHM_filer = FWHM_filer, flux_filter= [50,20000],
                                  nearyby_obj_filter=False, FWHM_sort=True)
            # data_process.plot_overview()
            PSFs = data_process.PSF_list
            # plt_many_fits(PSFs)
            # PSF_list = data_process.PSF_list
            PSF_pos_list = data_process.PSF_pos_list
            fov_PSF_pos_list = [PSF_pos_list[k]+np.array([ct_pos_j, ct_pos_i]) for k in range(len(PSF_pos_list))]
            all_PSF_pos_list_ = all_PSF_pos_list_ + fov_PSF_pos_list
    
    #%%
    cata_folder = 'Chandra_COSMOS_Catalog/'
    cata_file = 'chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits'
    hdul = pyfits.open(cata_folder+cata_file)
    table = hdul[1].data
    name = hdul[1].columns
    cata_list = []
    frame_flux = []
    target_pos_list = []
    for i in range(len(table)):
        RA, Dec = table[i][24], table[i][25]
        if RA != -99:
            pos = wcs.all_world2pix([[RA, Dec]], 1)[0]
            if pos[0]>0 and pos[1]>0 :
                try:
                    flux = img[int(pos[1]), int(pos[0]) ]  #!!! To be confirm if pos[1], pos[0]
                    if flux!=0: 
                        target_pos_list.append(pos)
                except:
                    continue
    target_pos_list = np.array(target_pos_list)
    remove_i = []
    for i in range(len(all_PSF_pos_list_)-1):
        if np.min(np.sqrt( np.sum( (all_PSF_pos_list_[i] - target_pos_list)**2,axis=1 ) )) < 10:
            remove_i.append(i)
            print('This PSF id{0} is repeated!'.format(i))
    
    
    #%%Remove the repeat
    all_PSF_pos_list = np.array(all_PSF_pos_list_)
    remove_i = []
    for i in range(len(all_PSF_pos_list)-1):
        if np.min(np.sqrt(np.sum( (all_PSF_pos_list[i] - all_PSF_pos_list[i+1:])**2,axis=1 ) )) < 10:
            remove_i.append(i)
            print('This PSF id{0} is repeated!'.format(i))
    
    #%%
    all_PSF_pos_list = [all_PSF_pos_list[i] for i in range(len(all_PSF_pos_list)) if i not in remove_i]
    
    for i,pos in enumerate(all_PSF_pos_list):
        pos[0], pos[1] = int(pos[0]), int(pos[1])
        # print(img[int(pos[1]),int(pos[0])])
        test_img = img[int(pos[1])-30:int(pos[1])+30, int(pos[0])-30:int(pos[0])+30]
        ct_pos = len(test_img)/2
        shift_pos = np.where(test_img == test_img.max())[0]-ct_pos, np.where(test_img == test_img.max())[1]-ct_pos, 
        pos[0] = pos[0]+shift_pos[1]
        pos[1] = pos[1]+shift_pos[0]
        
    #Clean up PSF:
    clean_up = True
    PSF_org_list = []
    PSF_clean_list = []
    PSF_RA_DEC_list = []
    if clean_up == True:
        # lines = np.loadtxt('target_info.txt', dtype='str')
        for pos in all_PSF_pos_list: 
            psf = img[int(pos[1])-50:int(pos[1])+50, int(pos[0])-50:int(pos[0])+50]
            PSF_org_list.append(psf)
            psf = psf_clean(psf,if_plot=False, nsigma=3, npixels=45, ratio_to_replace=0.005,
                            if_print_fluxratio=True)
            RA, Dec = wcs.all_pix2world([[pos[0], pos[1]]], 1)[0]
            PSF_clean_list.append(psf)
            PSF_RA_DEC_list.append([RA, Dec])
            # print("Before remove candidates")
        # plt_many_fits(PSF_org_list)
        plt_many_fits(PSF_clean_list)
        # pickle.dump([PSF_org_list, PSF_clean_list, all_PSF_pos_list, PSF_RA_DEC_list],
        #             open('material/'+filt+'_PSF_Library.pkl', 'wb'))

#%%Refine the PSF
# new_PSF_org_list, new_PSF_clean_list, new_all_PSF_pos_list, new_PSF_RA_DEC_list = [],[],[],[]
# from galight.tools.cutout_tools import cutout
# for i,pos in enumerate(all_PSF_pos_list):
#     image = cutout(image = img, center = pos, radius=120)
#     # plt_fits(image)
#     psf = psf_clean(image,if_plot=False, nsigma=3, npixels=45, ratio_to_replace=0.005,
#                     if_print_fluxratio=True)
#     plt_fits(psf)
#     ifsave = input('input any string to not save:\n')
#     if ifsave == '':
#         new_PSF_org_list.append(image)
#         new_PSF_clean_list.append(psf)
#         new_all_PSF_pos_list.append(all_PSF_pos_list[i])
#         new_PSF_RA_DEC_list.append(PSF_RA_DEC_list[i])
# pickle.dump([new_PSF_org_list, new_PSF_clean_list, new_all_PSF_pos_list, new_PSF_RA_DEC_list],
#             open('material/'+filt+'_PSF_Library_v2.pkl', 'wb'))
# plt_many_fits(new_PSF_clean_list)
#%%
print("After remove candidates")
PSF_lib_files = glob.glob('material/'+filt+'_PSF_Library.pkl')[0]
PSF_org_list, PSF_clean_list, all_PSF_pos_list, PSF_RA_DEC_list = pickle.load(open(PSF_lib_files,'rb'))
plt_many_fits(PSF_clean_list)