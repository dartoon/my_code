#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 10:45:59 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob
#RA, DEC information of the QSO:
QSO_RA, QSO_DEC = 35.2735048,  -4.6837584
#Calculate the zeropoint for HSC filters:
# file_header0 = fitsFile[0].header
# FLUXMAG0 = file_header0['FLUXMAG0']
zp = 27 # 2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import pickle
folders = glob.glob('20*-*/')
folders.sort()
import copy

#%%
# name = 'PSPS+sersic_fixpos_result'
# folders = ['2021-10-31']
# for folder in folders[:25]:
for folder in folders: 
    folder = folder[:-1]
    files = glob.glob(folder+'/*/*/*/*/*fits')
    for file in files:
        HSCband = file.split('calexp-')[1][:5]
        band = HSCband[-1]
        
        if glob.glob('{0}-band{1}.fits'.format(folder,band)) !=[]:
            print("already:", '{0}-band{1}'.format(folder,band))
            # continue
        else:
            print("cutout:", '{0}-band{1}'.format(folder,band))
            fitsFile = pyfits.open(file)
            #Load the fov image data:
            fov_image = fitsFile[1].data # check the back grounp
            #Derive the header informaion, might be used to obtain the pixel scale and the exposure time.
            header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
            #Derive the fov noise level map:
            err_data= fitsFile[3].data ** 0.5
            #Load the PSF data:
            PSF_files = glob.glob('psf*{0}*.fits'.format(HSCband))
            PSF = pyfits.getdata(PSF_files[0])
            
            data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [QSO_RA, QSO_DEC],
                                        pos_type = 'wcs', header = header,
                                      rm_bkglight = True, if_plot=False, zp = zp)
            
            #Generate the fitting materials
            data_process.generate_target_materials(radius=120, create_mask = False, nsigma=2.8,
                                                  exp_sz= 1.5, npixels = 40, if_plot=False)
            
            # data_process.apertures = [] #Assuming there is no host (i.e., as constant.) #!!!
    
            #Manually input the PSF:
            data_process.PSF_list = [PSF]
            print('{0}-band{1}'.format(folder,band))
            
            file_header = copy.deepcopy(fitsFile[1].header)
            file_header['CRPIX1'] = file_header['CRPIX1']-data_process.target_pos[0]+len(data_process.target_stamp)/2
            file_header['CRPIX2'] = file_header['CRPIX2']-data_process.target_pos[1]+len(data_process.target_stamp)/2
            
            pyfits.PrimaryHDU(data_process.target_stamp,header=file_header).writeto('{0}-band{1}.fits'.format(folder,band) ,overwrite=True)
            pyfits.PrimaryHDU(data_process.noise_map,header=file_header).writeto('{0}-band{1}_noise.fits'.format(folder,band) ,overwrite=True)
