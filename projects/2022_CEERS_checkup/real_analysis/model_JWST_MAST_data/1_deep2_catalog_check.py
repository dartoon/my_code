#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:17:34 2022

@author: Dartoon

Check the nature of the sources by comparing the deep 2 catalog
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

hdul = pyfits.open('../catalog_regions/zcat.deep2.dr4.uniq.fits')
data = hdul[1].data
cols = hdul[1].columns
names = cols.names
RAs, DECs = data['RA'], data['DEC']

#quick check if dis is can be calculated.
RA, Dec = 214.75522, 52.836795
dis_ = np.sqrt((RA-RAs)**2 + (Dec - DECs)**2)
dis = np.min(dis_)*3600  #arcsec
# if dis < 10:
#     idx = np.where(dis_ == np.min(dis_))[0][0]
#     print('\t',data['CLASS'][idx], data['ZBEST'][idx], data['ZERR'][idx], np.round(dis))

#%%
# filt = 'f356w'
# import glob
# folder = '/Users/Dartoon/Downloads/CEERS_JWST_data'
# filenames = glob.glob(folder+'/bkg_removed/*'+filt+'*.fits')
# from def_functions import target_in_fits
# targets_info = target_in_fits(filenames)
# for target in targets_info:
#     ID, RA, Dec, z_best = target[:4]
#     print('ID', ID, z_best[-1])
#     dis_ = np.sqrt((RA-RAs)**2 + (Dec - DECs)**2)
#     dis = np.min(dis_)*3600  #arcsec
#     # if dis < 10:
#     idx = np.where(dis_ == np.min(dis_))[0][0]
#     print('\t',data['CLASS'][idx], data['ZBEST'][idx], data['ZERR'][idx], np.round(dis))
        

#%%
hdul = pyfits.open('../catalog_regions/hlsp_candels_hst_wfc3_egs_v1_mass_cat.fits')
data = hdul[1].data
cols = hdul[1].columns
names = cols.names
RAs, DECs = data['RAdeg'], data['DECdeg']

# dis_ = np.sqrt((RA-RAs)**2 + (Dec - DECs)**2)
# dis = np.min(dis_)*3600  #arcsec
# if dis < 10:
#     idx = np.where(dis_ == np.min(dis_))[0][0]
#     print(data['Class_star'][idx], data['AGNflag'][idx], data['zbest'][idx])
filt = 'f356w'
import glob
folder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_CEERS/CEERS_JWST_data'
filenames = glob.glob(folder+'/bkg_removed/*'+filt+'*.fits')
from def_functions import target_in_fits
targets_info = target_in_fits(filenames)
#%%
for target in targets_info:
    ID, RA, Dec, z_best = target[:4]
    print('ID', ID, z_best[0], z_best[-1])
    dis_ = np.sqrt((RA-RAs)**2 + (Dec - DECs)**2)
    dis = np.min(dis_)*3600  #arcsec
    # if dis < 10:
    idx = np.where(dis_ == np.min(dis_))[0][0]
    print('\t', data['AGNflag'][idx],  round(dis,3), data['zbest'][idx], data['zspec'][idx])
        