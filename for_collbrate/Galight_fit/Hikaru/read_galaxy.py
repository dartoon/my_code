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

bands = 'GRIZY'
for i in [0]:
    for band in bands:
        idx = i
        filename = glob.glob('Radio_DOGs/{0}_*/*band-{1}.pkl'.format(idx, band))
        if filename == []:
            continue
        filename = filename[0]
        result = pickle.load(open(filename,'rb'))
        #%% Host data:
        data = result.fitting_specify_class.kwargs_data['image_data']    
        model = result.flux_2d_out['model']
        noise = result.fitting_specify_class.kwargs_data['noise_map']
        ID = filename.split('/')[1].split('_')[1]
        # band = filename.split('band-')[1][0]
        result.plot_final_qso_fit(target_ID = ID+'-'+ band)
        print(result.final_result_galaxy)
    
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
