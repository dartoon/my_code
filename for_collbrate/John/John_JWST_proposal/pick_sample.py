#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 16:36:47 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
# folder = './'
# folder = '/gpfs02/work/junyao.li/mock_quasar_simulation_z*/'
# names = 'sampling_gamma0.00_sigma0.30_run*.fits'
# filenames = glob.glob(folder+names)
# # print('filenames', filenames)
# filenames.sort()

# savefolder = './'
# i = 0
# for file in filenames:
#     # file = file.replace('=','\=')
#     filename = file.split('/')[-1]
#     print(i,file)
#     fitsfile = fits.open(file)
#     data = fitsfile[1].data
#     newdata = data[data['muv_quasar'] < -22]
#     hdu = pyfits.BinTableHDU(data=newdata)
#     hdu.writeto(savefolder+filename.split('.fits')[0]+'_selected.fits')
#     i = i+1
    
#%%
filenames = glob.glob('*selected.fits')
filenames.sort()

i = 0
for file in filenames:
    fitsfile = fits.open(file)
    if i == 0:
        data = fitsfile[1].data
        hdu = pyfits.BinTableHDU(data=data)
        hdu.writeto('Junyao_sim_selected.fits', overwrite=True)
    if i >=0 :
        fitsfile_tosave = fits.open('Junyao_sim_selected.fits')
        nrows1 = fitsfile[1].data.shape[0]
        nrows2 = fitsfile_tosave[1].data.shape[0]
        nrows = nrows1 + nrows2
        hdu = fits.BinTableHDU.from_columns(fitsfile_tosave[1].columns, nrows=nrows)
        for colname in fitsfile_tosave[1].columns.names:
            hdu.data[colname][nrows2:] = fitsfile[1].data[colname]
        hdu.writeto('Junyao_sim_selected.fits', overwrite=True)
    i = i+1
    print(i)
    # file = file.replace('=','\=')
