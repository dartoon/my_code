#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 16:50:19 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'/Users/Dartoon/Astro/Projects/Lens_Model_challenge/TDSLMC/simulating/material/real_source')
from source_info import  source_list #[0]: file name [1]: total size [2]: galfit R_e [3]:R_e/totalzise

for source_id in range(25):
# source_id = 4
    source_name=source_list[0][source_id]
    Re = source_list[2][source_id]
    print(source_name, source_id)
    hdu = pyfits.open('/Users/Dartoon/Astro/Projects/Lens_Model_challenge/TDSLMC/simulating/material/real_source/fix/{0}_fix.fits'.format(source_name))
    hd_gal_img = hdu[0].data
    hdu.close()
    print("source name:", source_name)
    from galight.tools.astro_tools import plt_fits
    plt_fits(hd_gal_img)
    
    # from skimage.transform import resize
    # new_image_0 = resize(hd_gal_img, [25,25])
    # plt_fits(new_image_0)
    
    # from scipy.ndimage import zoom
    # new_image_1 = zoom(hd_gal_img, 0.1)
    # plt_fits(new_image_1)
