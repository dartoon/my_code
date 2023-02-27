#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:02:07 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import glob
cata_list = pickle.load(open('material/cata_list.pkl','rb'))

name_list = ['cid_499', 'cid_509', 'cid_510', 'cid_531']
filts = ['F115W', 'F150W','F277W', 'F444W']
for idx in range(len(cata_list)):
    pos = cata_list[idx][3:5]
    target_name = cata_list[idx][-1]
    if target_name in name_list:
        print(idx, target_name)
        for filt in filts:
            files = glob.glob('fit_material/fit2_notrunyet_{1}_idx{0}*pkl'.format(idx,filt))
            if target_name == name_list[0]:
                for file in files:
                    psf_id = file.split('psf')[1].split('.pkl')[0]
                    fit_run = pickle.load(open(file,'rb'))
            #         psf = fit_run.fitting_specify_class.data_process_class.PSF_list[0]
            #         pyfits.PrimaryHDU(psf).writeto('Takumi_fitting_fits/PSFs/{0}_psf{1}.fits'.format(filt,psf_id),overwrite=True)
            # data_image = fit_run.fitting_specify_class.data_process_class.target_stamp
            # noise_image = fit_run.fitting_specify_class.data_process_class.noise_map
            # pyfits.PrimaryHDU(data_image).writeto('Takumi_fitting_fits/{1}_{0}_data.fits'.format(filt,target_name),overwrite=True)
            # pyfits.PrimaryHDU(noise_image).writeto('Takumi_fitting_fits/{1}_{0}_noise.fits'.format(filt,target_name),overwrite=True)
            #         # psf = 
            # # if '499' in target_name