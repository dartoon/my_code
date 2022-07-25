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

cata_file = '../catalog_regions/SDSS_DR16Q_v4.fits'
hdul = pyfits.open(cata_file)
table = hdul[1].data
name = hdul[1].columns
table = table[table['RA']>213]
table = table[table['RA']<216]
table = table[table['Dec']>52]
table = table[table['Dec']<54]

#%%
from astropy.wcs import WCS
import glob

folder = '/Users/Dartoon/Downloads/CEERS_JWST_data'
# files = '/*clear*/*f356w_i2d.fits'
files = '/*clear*/*f150w_i2d.fits'
data_file = glob.glob(folder+files)

targets = []
for i in range(len(data_file)):
    file = data_file[i]
    print(file)
    fitsFile = pyfits.open(file)
    fov_image = fitsFile[1].data # check the back grounp
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    wcs = WCS(header)
    for j in range(len(table)):
        RA, Dec,z = table[j][1], table[j][2],  table[j]['Z']
        pos = wcs.all_world2pix([[RA, Dec]], 1)[0]
        try:
            if pos[0]>0 and pos[1]>0 and fov_image[int(pos[1])-1, int(pos[0])-1 ] > 0 :
                if z >2:#  or z_photo>2 and z_spec ==-99.0:
                    print( table[j]['SDSS_NAME'], RA, Dec, fov_image[int(pos[1])-1, int(pos[0])-1 ], z )
        except:
            None


