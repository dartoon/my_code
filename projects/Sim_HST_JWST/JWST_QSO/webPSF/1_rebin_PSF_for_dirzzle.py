#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 16:11:45 2020

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
from astropy.cosmology import FlatLambdaCDM
import glob
import os

sys.path.insert(0, '../../share_tools')
#%%Set up data basic information
filt_l = ['F444W', 'F356W', 'F200W', 'F150W']
filt_id = 3
filt  = filt_l[filt_id]
rebin_folder_name = 'drizzle_PSF_'+ filt
os.mkdir(rebin_folder_name)
#%%rebin each PSF ID
for psf_id in range(9):
    psf_name = 'PSF_' + filt + '/PSF_id{0}.fits'.format(psf_id)
    print(psf_name)
    psf = pyfits.open(psf_name)
    psf_data = psf[0].data
    psf_data = psf_data[1:,1:]
    psf_data /= psf_data.sum()
    ##==============================================================================
    ## #Bin the image res. from high to low. 
    ##==============================================================================
    import rebin #From my share_tools
    factor=4
    pattern_x=[0,2,0,2,1,3,1,3]
    pattern_y=[0,0,2,2,3,3,1,1]      #from the info. given by observation
    ################Bin the PSF and save it################
    #exp_psf=rebin.expend_grid(psf_pixel_high_res)
    cut_fd=int((len(psf_data)-((int(len(psf_data)/8*2)-1)*4+3))/2)
    exp_psf_o=psf_data[1+cut_fd:-cut_fd,1+cut_fd:-cut_fd]+ 0  # To change it from 251 to 247.
    exp_psf=rebin.expend_grid(exp_psf_o)
    cut_len=int(round(len(exp_psf_o)/factor)*factor)
    cut_out_psf=np.zeros([len(pattern_x),cut_len,cut_len])
    image_bin_psf=np.zeros([len(pattern_x),int(cut_len/factor),int(cut_len/factor)])
    for i in range(len(pattern_x)):
        cut_out_psf[i]=exp_psf[pattern_x[i]:cut_len+pattern_x[i],pattern_y[i]:cut_len+pattern_y[i]]   #the size before bin
        image_bin_psf[i]=rebin.block(cut_out_psf[i],(int(cut_len/factor),int(cut_len/factor)),factor=factor)
        image_bin_psf[i] /= np.sum(image_bin_psf[i])  #unify the psf value
        pyfits.PrimaryHDU(image_bin_psf[i]).writeto(rebin_folder_name+'/non_drizzled_psf_id{0}-{1}.fits'.format(psf_id, i+1),overwrite=False)
#
