#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 16:31:40 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

idx = 1

from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

folder = 'gsf_id_{0}/'.format(idx)

steller_file = glob.glob(folder+'/gsf_spec_*.fits')[0]
hdul = pyfits.open(steller_file)
info_muv = hdul[1].header 

# steller_file = glob.glob(folder+'/SFH_*.fits')[0]
# hdul = pyfits.open(steller_file)
# info1 = hdul[0].header 


print(info_muv['MUV16'], info_muv['MUV50'], info_muv['MUV84'])