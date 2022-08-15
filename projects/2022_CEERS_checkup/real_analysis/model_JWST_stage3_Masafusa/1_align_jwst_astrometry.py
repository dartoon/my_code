#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 00:26:46 2022

@author: Dartoon

Load the CEERS HST image for quick check

Generate the jwst_shift_list.pkl.
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import warnings
warnings.filterwarnings("ignore")

from galight.data_process import DataProcess
import sys
sys.path.insert(0,'..')
from def_functions import RA_Dec_in_fit
from galight.tools.astro_tools import plt_fits
from galight.tools.measure_tools import detect_obj

#load HST
folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_HST_data/'
all_files= glob.glob(folder+'/egs_all_wfc3_ir_f160w_030mas_v1.9_drz.fits')  #For NIRCam
fitsFile = pyfits.open(all_files[0])
fov_image_HST = fitsFile[0].data # check the back grounp
header_HST = fitsFile[0].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
# plt_fits(fov_image[14000:14000+8000,19000:19000+8000])

#%%
#Load JWST
folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_JWST_Masafusa'
jwst_all_filenames = glob.glob(folder+'/bkg_removed/'+'*.fits')
f = open("material/target_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
result_folder = 'fit_result/'
lines = lines[1:]
# for line in enumerate(lines[2:]):
#%%
recording_list = []
show_print = True
# for idx in range(len(lines)):
for idx in range(1,2):
    line = lines[idx]
    target_id, RA, Dec, spec_z, photo_z = line.split(' ')
    RA, Dec, spec_z, photo_z = float(RA), float(Dec), float(spec_z), float(photo_z)
    # target_id = 'aegis_533'
    # RA, Dec = 214.87553, 52.866464
    data_process = DataProcess(fov_image = fov_image_HST, target_pos = [RA, Dec], pos_type = 'wcs', header = header_HST,
                              rm_bkglight = False, if_plot=False, zp = 27 )
    data_process.generate_target_materials(radius=100, create_mask = False, nsigma=2.8, if_select_obj=False,
                                          exp_sz= 1.2, npixels = 150, if_plot=False, cut_kernel = None)
    data_process.target_stamp = np.flip( data_process.target_stamp)
    apertures, segm_deblend, mask_apertures, tbl = detect_obj(data_process.target_stamp,
                                                              nsigma=2.8, npixels = 150)
    data_process.apertures = apertures
    rest_tbl = tbl
    res_x = rest_tbl['xcentroid'][rest_tbl['kron_flux'] == np.max(rest_tbl['kron_flux'])][0]
    res_y = rest_tbl['ycentroid'][rest_tbl['kron_flux'] == np.max(rest_tbl['kron_flux'])][0]
    if show_print==True:
        data_process.plot_aperture()
        print("centers:",round(res_x), round(res_y),'ap used', 
              rest_tbl['label'][rest_tbl['kron_flux'] == np.max(rest_tbl['kron_flux'])][0] )
    
    jwst_filenames = RA_Dec_in_fit(all_files=jwst_all_filenames, RA=float(RA), Dec=float(Dec))
    filters = [jwst_filenames[i].split('NIRCam')[1][2:7] for i in range(len(jwst_filenames))]
    jwst_filenames = [x for _,x in sorted(zip(filters,jwst_filenames))]
    filters.sort()
    
    shift_list_l = []
    shift_list_s = []
    shift_list = []
    for i, file in enumerate(jwst_filenames):
        _fitsFile = pyfits.open(file)
        fov_image = _fitsFile[1].data # check the back grounp
        header = _fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        # flux_mjsr = header['PHOTMJSR']
        if _fitsFile[0].header['CHANNEL'] == 'LONG':
            expsize = 1
        else:
            expsize = 2
        fov_noise_map = _fitsFile[2].data 
        data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
                                  rm_bkglight = False, if_plot=False, zp = 27, fov_noise_map=fov_noise_map )
        data_process.generate_target_materials(radius=100*expsize, create_mask = False, nsigma=2.8, if_select_obj=False,
                                              exp_sz= 1.2, npixels = 150*expsize, if_plot=False, cut_kernel = None)
        raw_tbl = data_process.tbl
        raw_tbl['kron_flux'] = np.nan_to_num(raw_tbl['kron_flux'])
        raw_x = raw_tbl['xcentroid'][raw_tbl['kron_flux'] == np.max(raw_tbl['kron_flux'])][0]
        raw_y = raw_tbl['ycentroid'][raw_tbl['kron_flux'] == np.max(raw_tbl['kron_flux'])][0]
        if show_print==True:
            data_process.plot_aperture(figsize=(4,3))
            print(filters[i], file.split('/')[-1])
            print(round(res_x*expsize-raw_x), round(res_y*expsize-raw_y), 'ap used',
                  raw_tbl['label'][raw_tbl['kron_flux'] == np.max(raw_tbl['kron_flux'])][0])
        shift = [round(res_x*expsize-raw_x), round(res_y*expsize-raw_y)]
        if _fitsFile[0].header['CHANNEL'] == 'LONG':
            shift_list_l.append(shift)
        else:
            shift_list_s.append(shift)
        shift_list.append(shift)
    warning = ''
    if shift_list_s != []:
        if np.max(np.max(shift_list_s,axis=0) - np.min(shift_list_s,axis=0))>15:
            warning = 'idx {0} {1}, shift SHORT not consistent.\n'.format(idx, target_id)
    if shift_list_l != []:
        if np.max(np.max(shift_list_l,axis=0) - np.min(shift_list_l,axis=0))>7:
            warning = 'idx {0} {1}, shift LONG not consistent.\n'.format(idx, target_id)
    if np.max(shift_list) > 30 or np.min(shift_list) < -30:
        warning = warning + ' shift too much.'
    recording_list.append([shift_list, [shift_list_s, shift_list_l], 
                          [jwst_filenames[i].split('/')[-1] for i in range(len(jwst_filenames))],
                          warning])
    print('idx {0} {1} finished...'.format(idx, target_id))
    
#%%
# for i in range(len(recording_list)):
#     if recording_list[i][-1] != '':
#         print(recording_list[i][-1])
# for i in range(len(recording_list)):
#     _file = 'NIRCam2_F200W_1_i2d_rmbkg.fits'
#     if _file in recording_list[i][2]:
#         idx = [idx for idx in range(len(recording_list[i][2]))if recording_list[i][2][idx] == _file][0]
#         print(_file, recording_list[i][0][idx])
ignore_id = [10, 21, 30, 31,  41, 46, 47, 52]
remove_id = [24, 55]
#%%
# import pickle
# # pickle.dump(recording_list, open('material/'+'jwst_shift_list.pkl', 'wb'))
# shift_list = pickle.load(open('material/jwst_shift_list.pkl','rb'))