#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 16:13:10 2022

@author: Dartoon

After 1_align_jwst_astromerty.py. 

Propose: Generate the data_process for JWST and HST together. Apertures will also be generated. 

Template: 1_test_model_allbands.py 
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

#%%
shift_list = pickle.load(open('material/jwst_shift_list.pkl','rb')) #For fov shift

#load HST
HST_folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_HST_data/'
HST_all_files= glob.glob(HST_folder+'/egs_all_wfc3_ir_f160w_030mas_v1.9_drz.fits')  #For NIRCam
HST_fitsFile = pyfits.open(HST_all_files[0])
fov_image_HST = HST_fitsFile[0].data # check the back grounp
header_HST = HST_fitsFile[0].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
# plt_fits(fov_image[14000:14000+8000,19000:19000+8000])

#%%
ignore_id = [10, 21, 30, 31,  41, 46, 47, 49, 52, 56]
remove_id = [24, 55]
#Load JWST
folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_JWST_Masafusa/bkg_removed/'
# jwst_all_filenames = glob.glob(folder+'/bkg_removed/'+'*.fits')
f = open("material/target_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
result_folder = 'fit_result/'
lines = lines[1:]
# for line in enumerate(lines[2:]):
# for idx in range(len(lines)):

cid = 57
# for idx in range(cid,cid+1):
for idx in range(cid,59):
    cut_kernel = None
    if idx in remove_id:
        continue
    line = lines[idx]
    target_id, RA, Dec, spec_z, photo_z = line.split(' ')
    RA, Dec, spec_z, photo_z = float(RA), float(Dec), float(spec_z), float(photo_z)
    files_list = shift_list[idx][2]
    # shift_list[idx][0][-1][0] = shift_list[idx][0][-1][0] - 10 #!!!
    # shift_list[idx][0][-1][1] = shift_list[idx][0][-1][1] + 5 #!!!
    # cut_kernel = 'nearest_obj_center'
    # filters = [files_list[i].split('NIRCam')[1][2:7] for i in range(len(files_list))]
    filters = []
    cut_RA, cut_Dec = RA, Dec
    data_process_list_l, data_process_list_s = [], []
    for i in range(len(files_list))[::-1]:  #Start with the reddest filter.
        cut_kernel = 'nearest_obj_center' #After pos correct then, do the nearest_obj_center
        # if i == 0:
        #     cut_kernel = None
        file = files_list[i]
        filt = file.split('NIRCam')[1][2:7]
        filters.append(filt)
        print("Loading...,", 'idx', idx, file)
        _fitsFile = pyfits.open(folder+file)
        fov_image = _fitsFile[1].data # check the back grounp
        header = _fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        # flux_mjsr = header['PHOTMJSR']
        flux_mjsr = header['PHOTMJSR']
        pixscale = read_pixel_scale(header)
        zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
        wht = _fitsFile[4].data # The WHT map
        exp = _fitsFile[0].header['EFFEXPTM']
        if _fitsFile[0].header['CHANNEL'] == 'LONG':
            gain_value = 2
            expsize = 1
            exppix = 1
        else:
            gain_value = 1.8
            expsize = 1.4
            exppix = 2
        exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
        # fov_noise_map = _fitsFile[2].data 
        print("Processing data...")
        data_process = DataProcess(fov_image = fov_image, target_pos = [cut_RA, cut_Dec], #The final cut center is dete
                                   pos_type = 'wcs', header = header,rm_bkglight = False, 
                                   if_plot=False, zp = zp, exptime= exp_map, 
                                   fov_noise_map = None)
        if idx not in ignore_id:
            data_process.target_pos = data_process.target_pos - np.array(shift_list[idx][0][i]) + np.array(shift_list[idx][0][-1]) * exppix
        #estimate local bkg and remove:
        data_process.generate_target_materials(radius=120 * expsize, skip = True,
                                                cut_kernel = cut_kernel, 
                                                if_plot=False, npixels = 300 * expsize)
        cut_kernel = None
        if idx in ignore_id:
            wcs = WCS(header)
            cut_RA, cut_Dec = wcs.all_pix2world([data_process.target_pos], 1)[0] #re-define RA, Dec
        bkg_std = data_process.bkg_std
        fov_cutout = cutout(image=data_process.fov_image, center= data_process.target_pos, radius=200 * expsize)
        
        bkglight = measure_bkg(fov_cutout, if_plot=False) # Remove bkg light
        
        data_process.generate_target_materials(radius=45 * expsize, create_mask = False, nsigma=2.0, 
                                                cut_kernel = None, if_select_obj=False,
                                              exp_sz= 1.2, npixels = 100 * expsize, if_plot=False, bkg_std= bkg_std)
        data_process.noise_map[data_process.noise_map == 0] = data_process.noise_map.max()
        if np.sum(data_process.target_stamp ==0) >5:
            data_process.target_mask = data_process.target_stamp != 0
            data_process.noise_map = np.nan_to_num(data_process.noise_map, nan=1000)
        ct = int((len(bkglight) - len(data_process.target_stamp ))/2)
        data_process.target_stamp = data_process.target_stamp - bkglight[ct:-ct, ct:-ct]

        del data_process.fov_image
        del data_process.exptime
        print('idx', idx, target_id, filt, 'apertures', len(data_process.apertures) )
        data_process.filt = filt
        data_process.file = file
        # data_process.plot_aperture()
        data_process.cut_RA, data_process.cut_Dec =  cut_RA, cut_Dec
        if _fitsFile[0].header['CHANNEL'] == 'LONG':
            data_process_list_l.append(data_process)
        if _fitsFile[0].header['CHANNEL'] == 'SHORT':
            data_process_list_s.append(data_process)
    com_aper_s, com_aper_l = [], []
    if data_process_list_l != []:
        com_aper_l = common_data_class_aperture(data_process_list_l, l_idx=0, return_idx=0)
    if data_process_list_s != []:
        com_aper_s = common_data_class_aperture(data_process_list_s, l_idx=0, return_idx=0)
    # com_aper_s = data_process_list_s[0].apertures
    # del com_aper_s[1]
    # del com_aper_s[3]
    # com_aper_l = data_process_list_l[0].apertures
    # com_aper_s = data_process_list_s[0].apertures
    
    print("Common aperture:")
    for i in range(len(data_process_list_l)):
        print(idx, target_id, data_process_list_l[i].filt,':')
        plot_data_apertures_point(data_process_list_l[i].target_stamp * data_process_list_l[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                                  com_aper_l, figsize=(4,3))
    for i in range(len(data_process_list_s)):
        print(idx, target_id, data_process_list_s[i].filt,':')
        plot_data_apertures_point(data_process_list_s[i].target_stamp * data_process_list_s[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                                  com_aper_s, figsize=(4,3))
    print("Above are for the", 'idx:', idx, target_id, 'filts:', filters)
    hold = input('Hold ... OK?\n')
    if hold == '!':
        break
    pickle.dump([[data_process_list_l, com_aper_l], [data_process_list_s, com_aper_s] ], open('material/'+'data_process+apertures_{0}.pkl'.format(idx), 'wb'))
    
    import shutil
    shutil.copyfile('1_generate_data_process.py', 'backup_files/1_generate_data_process_idx{0}.py'.format(idx))
