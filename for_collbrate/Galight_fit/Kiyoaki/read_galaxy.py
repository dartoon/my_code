#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 21:51:36 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle
import glob

for i in range(20):
    idx = i
    filename = glob.glob('{0}_*/*pkl'.format(idx))
    if filename == []:
        continue
    filename = filename[0]
    result = pickle.load(open(filename,'rb'))
    #%% Host data:
    data = result.fitting_specify_class.kwargs_data['image_data']    
    model = result.flux_2d_out['model']
    noise = result.fitting_specify_class.kwargs_data['noise_map']
    ID = filename.split('/')[0].split('_')[1]
    result.plot_final_galaxy_fit(target_ID = ID)

    # %%
    # from matplotlib.colors import LogNorm
    # norm = LogNorm()
    # plt.imshow(data, origin='lower', norm = norm) 
    # plt.colorbar()
    # plt.show() 
    
    # plt.imshow(model, origin='lower', norm = norm) 
    # plt.colorbar()
    # plt.show() 
    
    # plt.imshow((data-model)/noise, origin='lower', norm = norm) 
    # plt.colorbar()
    # plt.show() 
