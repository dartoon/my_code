#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 10:24:56 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

files = ['2016-02-09/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-12-28/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-12-31/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-09-29/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2020-01-20/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-10-27/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2020-01-26/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-11-02/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2020-02-22/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2020-02-24/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-11-27/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2021-09-01/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-12-01/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-12-02/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2021-10-29/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-12-04/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-12-22/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits']

path = '/gpfs02/work/yasuda/transient/sxds2019/rerun/sxds2019/rerun/'
import shutil
import os
for file in files:
    folder = file.split('cal')[0]
    print(folder)
    os.makedirs(os.path.dirname(folder))
    src_fpath = path+file
    dest_fpath = os.path.dirname(folder)
    shutil.copy(src_fpath, dest_fpath)
