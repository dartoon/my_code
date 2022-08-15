#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:53:44 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from galight.tools.astro_tools import plt_fits

f = open('../catalog_regions/AEGIS_data_140612.csv',"r") ##This RA DEC of the optical counterparts is used to get hte reg file
string = f.read()
AEGIS_2014 = string.split('\n')   # Split in to \n

f = open('../catalog_regions/AEGIS-XD_redshift_catalog.txt',"r") ##This RA DEC of the optical counterparts is used to get hte reg file
string = f.read()
AEGIS_redshift = string.split('\n')   # Split in to \n

#%%
from astropy.wcs import WCS
import glob

def return_z(target_id):
    line = [AEGIS_redshift[i] for i in range(len(AEGIS_redshift)) if target_id in AEGIS_redshift[i]]
    if len(line)>1:
        print("find two same ID in redshift.")#, line)
    s = line[0]
    info = s.split(' ')
    info = [info_ for info_ in info if  info_ !='']
    if info[4] != '0' and float(info[3])>=0:
        z_spec = float(info[2])
    else:
        z_spec = -99.
    z_photo = float(info[6])
    return z_spec, z_photo
folder = '/Users/Dartoon/Downloads/CEERS_JWST_data'
# files = '/*clear*/*f356w_i2d.fits'
files = '/*clear*/*f150w_i2d.fits'
data_file = glob.glob(folder+files)

targets = []
for i in range(len(data_file)):
    file = data_file[i]
    fitsFile = pyfits.open(file)
    fov_image = fitsFile[1].data # check the back grounp
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    # plt_fits(fov_image)
    wcs = WCS(header)
    # print('search for file', file)
    for line in AEGIS_2014[17:-1]:
        info = line.split(' ')
        info = [info_ for info_ in info if  info_ !=''] 
        target_id, RA, Dec = info[0], float(info[3]), float(info[4])
        
        if RA != -99:
            pos = wcs.all_world2pix([[RA, Dec]], 1)[0]
        try:
            if pos[0]>0 and pos[1]>0 and fov_image[int(pos[1])-1, int(pos[0])-1 ] > 0 :
                z_spec, z_photo = return_z(target_id)
                if z_spec >2:# or z_photo>2 and z_spec ==-99.0:
                    print(file)
                    print(target_id, RA, Dec, 'flux:',fov_image[int(pos[1])-1, int(pos[0])-1 ], 'redshift:', z_spec, z_photo)
                    # targets.append([target_id, RA, Dec, 'flux:',fov_image[int(pos[1])-1, int(pos[0])-1 ], 'redshift:', z_spec, z_photo])
        except:
            None
# #The same data in these fits (i.e., t21):
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-c1002_t021_nircam_clear-f356w/jw01345-c1002_t021_nircam_clear-f356w_i2d.fits
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o001_t021_nircam_clear-f356w/jw01345-o001_t021_nircam_clear-f356w_i2d.fits
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-c1001_t021_nircam_clear-f356w/jw01345-c1001_t021_nircam_clear-f356w_i2d.fits
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-c1000_t021_nircam_clear-f356w/jw01345-c1000_t021_nircam_clear-f356w_i2d.fits
        
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-c1001_t021_nircam_clear-f150w/jw01345-c1001_t021_nircam_clear-f150w_i2d.fits
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-c1002_t021_nircam_clear-f150w/jw01345-c1002_t021_nircam_clear-f150w_i2d.fits
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o001_t021_nircam_clear-f150w/jw01345-o001_t021_nircam_clear-f150w_i2d.fits
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-c1000_t021_nircam_clear-f150w/jw01345-c1000_t021_nircam_clear-f150w_i2d.fits

# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o003_t023_nircam_clear-f356w/jw01345-o003_t023_nircam_clear-f356w_i2d.fits
# aegis_463 214.7768 52.825876 flux: 1.0515807 redshift: 2.274 2.24
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o003_t023_nircam_clear-f356w/jw01345-o003_t023_nircam_clear-f356w_i2d.fits
# aegis_482 214.75522 52.836795 flux: 0.9378855 redshift: 3.465 3.345
# /Users/Dartoon/Downloads/CEERS_JWST_data/jw01345-o004_t024_nircam_clear-f356w/jw01345-o004_t024_nircam_clear-f356w_i2d.fits
# aegis_477 214.87073 52.833117 flux: 0.5838029 redshift: 2.317 2.136
