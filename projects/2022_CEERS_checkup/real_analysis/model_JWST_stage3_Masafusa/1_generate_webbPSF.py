#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 21:45:00 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt


import webbpsf
oversample = 2
nc = webbpsf.NIRCam()
nc.detector = 'NRCA5'
filt_list = ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']
for filt in filt_list:
    nc.filter = filt
    POS = [1000, 1000]
    if int(filt[1:2])>=21:
        nc.detector = 'NRCB5'
    else:
        nc.detector = 'NRCB1'
    nc.detector_position = POS
    psf_webb_fits = nc.calc_psf(oversample=oversample)
    psf_webb = psf_webb_fits[0].data
    pyfits.PrimaryHDU(psf_webb).writeto('webbPSFs/PSF_{0}.fits'.format(filt),overwrite=True)
    print('filt', filt, "finished")