#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 14:10:06 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

filt = ['F115W', 'F150W','F277W', 'F444W'][1]
filefolder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/'
filename = 'mosaic_nircam_f{0}w_COSMOS-Web_30mas_v0_1_i2d.fits'.format(filt[1:-1])
# filename = '1727_cosmos_mosaic_miri_exptime_scale1.0.fits'
fitsFile = pyfits.open(filefolder+filename)
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
img = fitsFile[1].data #
print(img.shape)
# print(img[14000,20000])
# from galight.tools.astro_tools import plt_fits
# plt_fits(img)

#%%
from astropy.wcs import WCS
wcs = WCS(header)
cata_folder = 'Chandra_COSMOS_Catalog/'
cata_file = 'chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits'
hdul = pyfits.open(cata_folder+cata_file)
table = hdul[1].data
name = hdul[1].columns
cata_list = []
frame_flux = []
from galight.tools.astro_tools import plt_fits
for i in range(len(table)):
    RA, Dec = table[i][24], table[i][25]
    if RA != -99:
        pos = wcs.all_world2pix([[RA, Dec]], 1)[0]
        if pos[0]>0 and pos[1]>0 :
            try:
                flux = img[int(pos[1]), int(pos[0]) ]  #!!! To be confirm if pos[1], pos[0]
                if flux!=0: 
                    print(pos, flux, table[i][38], table[i][39])
                    cata_list.append([i, RA, Dec, pos[0], pos[1], table[i][38], table[i][39], table[i][0]]) #myid, RA, Dec, pos, redshift, name 
                    plt_fits(img[int(pos[1])-100:int(pos[1])+100, int(pos[0])-100:int(pos[0])+100])
                    # flux = np.sum(img[int(pos[1])-40:int(pos[1])+40, int(pos[0])-40:int(pos[0])+40])
                    print(flux)
                    frame_flux.append(flux)
            except:
                continue
    # cata_list = np.array(cata_list)
#%%
# import pickle
# pickle.dump(cata_list, open('material/'+'cata_list.pkl', 'wb'))