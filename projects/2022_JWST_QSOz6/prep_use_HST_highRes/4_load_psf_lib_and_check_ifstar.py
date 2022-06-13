#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 17:22:06 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy
import pickle
filt = 'f356w'

# psfs, FWHMs = pickle.load(open('f150w_psfs.pkl','rb'))
psfs, FWHMs = pickle.load(open(filt+'_psfs.pkl','rb'))

# from galight.tools.astro_tools import plt_many_fits
# for i in range(len(psfs)):
#     # plt_many_fits(psfs[i], FWHMs[i], 'FWHM')
#     plt_many_fits(psfs[i], [np.sum(psfs[i][j]) for j in range(len(psfs[i]))], 'Flux')
# #     if len(psfs[i]) != len(FWHMs[i]):   
# #           print(i, "Shape not right for!")

folder = 'JWST_CEERS/'
file = 'ceers5_{filt}_i2d.fits'.format(filt=filt)
im = pyfits.open(folder+file)
data = im[1].data
header = im[1].header
# for i in range(len(psfs[0])):
#     psf = psfs[0][i]
#     print(np.where(psf.max() == psf))

PS_list = np.loadtxt(folder+'ptsrc.cat')
PS_list = PS_list[:,1:3]

PS_pix_list = []
from astropy.wcs import WCS
for i in range(len(PS_list)):
    Ra, Dec = PS_list[i]
    wcs = WCS(header)
    target_pos = wcs.all_world2pix([[Ra, Dec]], 1)[0]
    PS_pix_list.append(target_pos)
    # print(target_pos, Ra, Dec)
PS_pix_list = np.array(PS_pix_list)
    
wcs = WCS(header)
for ii, psfs_ in enumerate(psfs[:]):
    for jj, psf in enumerate(psfs_[:]): 
        data_ = copy.deepcopy(data)
        x, y = np.where(psf.max() == psf)
        
        s_ref = psf[x,y]/psf[x-1,y]
        for i in range(10):
            find_pos = np.where( np.min(abs(data_-psf[x,y]) ) ==  abs(data_-psf[x,y]) )
            x_, y_ = find_pos[0][0], find_pos[1][0]
            f_ref = data_[x_,y_]/data[x_-1,y_]
            if abs(s_ref - f_ref) < 5.e-4:
                # print( data_[x_, y_], psf.max(), f_ref, s_ref[0] , i) 
                break
            else:
                data_[x_,y_] = -99
        print(ii, jj)
        if np.sqrt(np.min((np.array(y_, x_) - PS_pix_list)**2)) <5:
            print("Point source", np.sqrt(np.min((np.array(y_, x_) - PS_pix_list)**2)))
        else:
            print("Not in catalog", np.sqrt(np.min((np.array(y_, x_) - PS_pix_list)**2)))
        
        
            