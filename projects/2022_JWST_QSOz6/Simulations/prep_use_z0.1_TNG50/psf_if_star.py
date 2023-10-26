#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 11:31:22 2022

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import copy

def if_star(psf,filt = 'f150w'):
    folder = '../prep_use_HST_highRes/JWST_CEERS/'
    file = 'ceers5_{filt}_i2d.fits'.format(filt=filt)
    im = pyfits.open(folder+file)
    data = im[1].data
    header = im[1].header
    wcs = WCS(header)
    data_ = copy.deepcopy(data)
    x, y = np.where(psf.max() == psf)
    
    PS_list = np.loadtxt(folder+'ptsrc.cat')
    PS_list = PS_list[:,1:3]
    
    PS_pix_list = []
    for i in range(len(PS_list)):
        Ra, Dec = PS_list[i]
        wcs = WCS(header)
        target_pos = wcs.all_world2pix([[Ra, Dec]], 1)[0]
        PS_pix_list.append(target_pos)
        # print(target_pos, Ra, Dec)
    PS_pix_list = np.array(PS_pix_list)
    
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
    dis_pix = np.sqrt(np.min((np.array(y_, x_) - PS_pix_list)**2))
    if dis_pix <5:
        print("Point source; dis pix ", round(dis_pix,2))
        return True
    else:
        print("Not in catalog; dis pix ", round(dis_pix,2))
        return False