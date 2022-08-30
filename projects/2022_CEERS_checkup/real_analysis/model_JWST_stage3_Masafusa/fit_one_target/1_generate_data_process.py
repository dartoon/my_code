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
import sys
sys.path.insert(0,'../..')
from def_functions import RA_Dec_in_fit

#%%
target_id, RA, Dec = 'SDSS_0', 214.8234583, 52.83030556
folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_JWST_Masafusa'
filenames = glob.glob(folder+'/bkg_removed/'+'*.fits')
files_list = RA_Dec_in_fit(all_files=filenames, RA=float(RA), Dec=float(Dec))
filters = [files_list[i].split('NIRCam')[1][2:7] for i in range(len(files_list))]
files = [x for _,x in sorted(zip(filters,files_list))]
# #Load JWST
# folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_JWST_Masafusa/bkg_removed/'

cut_kernel = None
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
    _fitsFile = pyfits.open(file)
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
    #estimate local bkg and remove:
    data_process.generate_target_materials(radius=55 * expsize, skip = True,
                                            cut_kernel = cut_kernel, 
                                            if_plot=False, npixels = 30 * expsize)
    print(data_process.target_pos)
    cut_kernel = None
    
    wcs = WCS(header)
    cut_RA, cut_Dec = wcs.all_pix2world([data_process.target_pos], 1)[0] #re-define RA, Dec
    bkg_std = data_process.bkg_std
    fov_cutout = cutout(image=data_process.fov_image, center= data_process.target_pos, radius=200 * expsize)
    
    bkglight = measure_bkg(fov_cutout, if_plot=False) # Remove bkg light
    from galight.tools.astro_tools import plt_fits
    plt_fits(fov_cutout)
    data_process.generate_target_materials(radius=30 * expsize, create_mask = False, nsigma=1.2, 
                                            cut_kernel = None, if_select_obj=False,
                                          exp_sz= 1.2, npixels = 40 * expsize, if_plot=False, bkg_std= bkg_std)
    data_process.noise_map[data_process.noise_map == 0] = data_process.noise_map.max()
    if np.sum(data_process.target_stamp ==0) >5:
        data_process.target_mask = data_process.target_stamp != 0
        data_process.noise_map = np.nan_to_num(data_process.noise_map, nan=1000)
    ct = int((len(bkglight) - len(data_process.target_stamp ))/2)
    data_process.target_stamp = data_process.target_stamp - bkglight[ct:-ct, ct:-ct]

    del data_process.fov_image
    del data_process.exptime
    print('target_id', target_id, filt, 'apertures', len(data_process.apertures) )
    data_process.filt = filt
    data_process.file = file
    data_process.plot_aperture()
    data_process.cut_RA, data_process.cut_Dec =  cut_RA, cut_Dec
    if _fitsFile[0].header['CHANNEL'] == 'LONG':
        data_process_list_l.append(data_process)
    if _fitsFile[0].header['CHANNEL'] == 'SHORT':
        data_process_list_s.append(data_process)
    # hold = input('OK?')

#%%    
com_aper_s, com_aper_l = [], []
if data_process_list_l != []:
    com_aper_l = common_data_class_aperture(data_process_list_l, l_idx=0, return_idx=0)
if data_process_list_s != []:
    com_aper_s = common_data_class_aperture(data_process_list_s, l_idx=0, return_idx=0)

print("Common aperture:")
for i in range(len(data_process_list_l)):
    print(target_id, data_process_list_l[i].filt,':')
    plot_data_apertures_point(data_process_list_l[i].target_stamp * data_process_list_l[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                              com_aper_l, figsize=(4,3))
for i in range(len(data_process_list_s)):
    print(target_id, data_process_list_s[i].filt,':')
    plot_data_apertures_point(data_process_list_s[i].target_stamp * data_process_list_s[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                              com_aper_s, figsize=(4,3))
print(target_id, 'filts:', filters)
hold = input('Hold ... OK?\n')
pickle.dump([[data_process_list_l, com_aper_l], [data_process_list_s, com_aper_s] ], open('material/'+'data_process+apertures_{0}.pkl'.format(target_id), 'wb'))

import shutil
shutil.copyfile('1_generate_data_process.py', 'backup_files/1_generate_data_process_{0}.py'.format(target_id))
