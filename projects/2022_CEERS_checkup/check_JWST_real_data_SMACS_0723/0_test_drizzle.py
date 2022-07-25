#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 13:25:27 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

raw_files = glob.glob('./F356W/*cal.fits')
#%% Remove bkg for each dither
for file in raw_files:
    fitsFile = pyfits.open(file, mode='update')
    fov_image = fitsFile[1].data # check the back grounp
    print(fov_image.max())
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    #
    wht = fitsFile[2].data # The WHT map
    import galight.tools.astro_tools as astro_tools
    exp =   header['XPOSURE'] #Read the exposure time 
    mean_wht = exp #* (0.0642/0.135)**2
    exp_map = exp * wht/mean_wht
    
    #Start to use galight
    from galight.data_process import DataProcess
    data_process = DataProcess(fov_image = fov_image, target_pos = [0., 0.], pos_type = 'pixel', header = header,
                              rm_bkglight = True, exptime = exp_map, if_plot=True, zp = 27.0)
    print(data_process.fov_image.max())
    fitsFile[1].data = data_process.fov_image
    fitsFile.flush()
    
# print("total drz:", len(files))
        