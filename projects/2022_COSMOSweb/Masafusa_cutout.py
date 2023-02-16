#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 14:57:25 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

folder = 'otherfiles/'
cata_file = 'COSMOSWeb_DEC_v0.23.5mf_prep2_z6can.fits'
# cata_file = 'COSMOSWeb_DEC_v0.23.5mf_prep2_cosmos2020.fits'
hdul = pyfits.open(folder+cata_file)
table = hdul[1].data
name = hdul[1].columns
name_list = [table[i][0] for i in range(len(table)) ]
RA_list = [table[i][1] for i in range(len(table)) ]
Dec_list = [table[i][2] for i in range(len(table)) ]
label_id = [i for i in range(len(table))]

#%%
filt_i = 3
filt = ['F115W', 'F150W','F277W', 'F444W', 'F770W'][filt_i]
filefolder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/'
if filt != 'F770W':
    filename = 'mosaic_nircam_f{0}w_COSMOS-Web_30mas_v0_1_i2d.fits'.format(filt[1:-1])
else:
    filename = 'mosaic_miri_f770w_COSMOS-Web_60mas_v0_1_i2d.fits'
fitsFile = pyfits.open(filefolder+filename)
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
img = fitsFile[1].data #
print(img.shape)
from astropy.wcs import WCS
wcs = WCS(header)

cut_p = 5
zp = 28
# pos = wcs.all_world2pix([[RA, Dec]], 1)[0]
from galight.data_process import DataProcess
from galight.tools.astro_tools import plt_fits, plt_many_fits
image_i = []
for i in range(len(RA_list)):
    name = name_list[i]
    RA, Dec = RA_list[i],Dec_list[i]
    data_process = DataProcess(fov_image = img, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
                              rm_bkglight = False, if_plot=False, zp = zp, fov_noise_map = img)
    try:
        data_process.generate_target_materials(radius=None, create_mask = False, nsigma=1, if_select_obj=False,
                                              exp_sz= 1.2, npixels = 5, if_plot=False)
        # plt_fits(data_process.target_stamp)
        image_i.append(data_process.target_stamp)
    except:
        image_i.append(np.ones([50,50]))
    # print(i)
#%%
plt_many_fits(image_i, labels = name_list, texts = label_id, prop = 'id', norm=None)