#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 17:16:06 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle, glob
import warnings
warnings.filterwarnings("ignore")
from matplotlib.colors import LogNorm
from galight.tools.astro_tools import read_pixel_scale

cata_list = pickle.load(open('../material/cata_list.pkl','rb'))
check_name= 'cid_473'  #29
# check_name= 'cid_1210' #8
# check_name= 'cid_1245' #10
check_id = [i for i in range(len(cata_list)) if cata_list[i][-1] == check_name]
print(cata_list[check_id[0]])

from galight.data_process import DataProcess
from scipy.ndimage import zoom
filts = ['F115W', 'F150W','F277W', 'F444W']
filefolder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/'
size = 90  #For LW it is 56*4 and for SW it is 56*2
sed_2d_info  = []
image_list = [None] * len(filts)
zp_dict = {}
for idx in check_id:
    for i, filt in enumerate(filts):
        filename = 'mosaic_nircam_f{0}w_COSMOS-Web_30mas_v0_1_i2d.fits'.format(filt[1:-1])
        fitsFile = pyfits.open(filefolder+filename)
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        img = fitsFile[1].data #
        pos = cata_list[idx][3:5]
        flux_mjsr = header['PHOTMJSR']
        pixscale = read_pixel_scale(header)
        zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
        data_process = DataProcess(fov_image = img, target_pos = pos, pos_type = 'pixel', header = header,
                                  rm_bkglight = False, if_plot=False, zp = 27, fov_noise_map = img)
        data_process.generate_target_materials(radius=size,create_mask = False, nsigma=2.8, if_select_obj=False,
                                              exp_sz= 1.2, npixels = 100, if_plot=False)

        img_show = zoom(data_process.target_stamp, 0.5)
        img_show = img_show/np.sum(img_show) * np.sum(data_process.target_stamp)
        image_list[i] = img_show
        zp_dict[filt] = zp
        print(filt,'finish', img_show.shape)
images = []
# zp_list = []
for i in [-1,-2,-4]:  #['F356W', 'F200W', 'F115W', 'F150W', 'F277W', 'F410M', 'F444W']
    # zp_list.append(zp_dict[filts[i]])
    images.append(image_list[i])
from galight.tools.astro_tools import plt_fits_color, plt_fits
# images = [images[i] * 10 ** (-0.4*(zp_list[i]-zp_list[0])) for i in range(len(images)) ]
for i in range(len(images)):
    plt_fits(images[i], vmin=0.001, vmax=2.5)
plt_fits_color(images, Q=7, stretch=0.3)


#%%        
target_id, z = cata_list[check_id[0]][-1], cata_list[check_id[0]][-2]
sed_image = np.zeros_like(image_list[0])
sed_2d_info = []
count = 1
for i in range(len(sed_image[0])):
    for j in range(len(sed_image[1])):
        mag_result = {}
        for k in range(len(filts)):
            filt = filts[k]
            flux = image_list[k][i,j]
            if flux>0:
                mag = -2.5*np.log10(flux) + zp_dict[filt]
                mag_result[filt] = mag
        sed_2d_info.append([i, j, mag_result])
import pickle      
pickle.dump(sed_2d_info, open('2d_filts_mag_bin2_{0}.pkl'.format(check_name), 'wb'))


