#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 10:24:56 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

files = ['2016-01-12/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-12-28/deepCoadd/HSC-G/8524/6,7/calexp-HSC-G-8524-6,7.fits', 
'2016-02-09/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-12-28/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-09-27/deepCoadd/HSC-G/8524/6,7/calexp-HSC-G-8524-6,7.fits', 
'2019-12-31/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-09-28/deepCoadd/HSC-G/8524/6,7/calexp-HSC-G-8524-6,7.fits', 
'2019-12-31/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-09-29/deepCoadd/HSC-G/8524/6,7/calexp-HSC-G-8524-6,7.fits', 
'2020-01-01/deepCoadd/HSC-Y/8524/6,7/calexp-HSC-Y-8524-6,7.fits', 
'2019-09-29/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2020-01-19/deepCoadd/HSC-G/8524/6,7/calexp-HSC-G-8524-6,7.fits', 
'2019-09-29/deepCoadd/HSC-R2/8524/6,7/calexp-HSC-R2-8524-6,7.fits', 
'2020-01-19/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-09-29/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2020-01-20/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-10-03/deepCoadd/HSC-Y/8524/6,7/calexp-HSC-Y-8524-6,7.fits', 
'2020-01-20/deepCoadd/HSC-R2/8524/6,7/calexp-HSC-R2-8524-6,7.fits', 
'2019-10-04/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2020-01-21/deepCoadd/HSC-Y/8524/6,7/calexp-HSC-Y-8524-6,7.fits', 
'2019-10-23/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2020-01-22/deepCoadd/HSC-G/8524/6,7/calexp-HSC-G-8524-6,7.fits', 
'2019-10-26/deepCoadd/HSC-R2/8524/6,7/calexp-HSC-R2-8524-6,7.fits', 
'2020-01-22/deepCoadd/HSC-R2/8524/6,7/calexp-HSC-R2-8524-6,7.fits', 
'2019-10-27/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2020-01-26/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-10-27/deepCoadd/HSC-R2/8524/6,7/calexp-HSC-R2-8524-6,7.fits', 
'2020-01-26/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-10-31/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2020-01-28/deepCoadd/HSC-Y/8524/6,7/calexp-HSC-Y-8524-6,7.fits', 
'2019-11-01/deepCoadd/HSC-G/8524/6,7/calexp-HSC-G-8524-6,7.fits', 
'2020-02-19/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-11-02/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2020-02-21/deepCoadd/HSC-Y/8524/6,7/calexp-HSC-Y-8524-6,7.fits', 
'2019-11-05/deepCoadd/HSC-Y/8524/6,7/calexp-HSC-Y-8524-6,7.fits', 
'2020-02-22/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-11-23/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2020-02-24/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-11-27/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2020-02-25/deepCoadd/HSC-R2/8524/6,7/calexp-HSC-R2-8524-6,7.fits', 
'2019-11-27/deepCoadd/HSC-R2/8524/6,7/calexp-HSC-R2-8524-6,7.fits', 
'2020-02-26/deepCoadd/HSC-G/8524/6,7/calexp-HSC-G-8524-6,7.fits', 
'2019-11-27/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2021-08-31/deepCoadd/HSC-G/8524/6,7/calexp-HSC-G-8524-6,7.fits', 
'2019-11-28/deepCoadd/HSC-G/8524/6,7/calexp-HSC-G-8524-6,7.fits', 
'2021-09-01/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-11-29/deepCoadd/HSC-R2/8524/6,7/calexp-HSC-R2-8524-6,7.fits', 
'2021-10-10/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-12-01/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2021-10-13/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-12-02/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2021-10-28/deepCoadd/HSC-R2/8524/6,7/calexp-HSC-R2-8524-6,7.fits', 
'2019-12-02/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2021-10-29/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2019-12-03/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2021-10-29/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-12-04/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2021-10-30/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-12-04/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2021-10-31/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-12-20/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2021-11-01/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-12-22/deepCoadd/HSC-I2/8524/6,7/calexp-HSC-I2-8524-6,7.fits', 
'2021-11-10/deepCoadd/HSC-Z/8524/6,7/calexp-HSC-Z-8524-6,7.fits', 
'2019-12-22/deepCoadd/HSC-R2/8524/6,7/calexp-HSC-R2-8524-6,7.fits', 
'2021-12-10/deepCoadd/HSC-G/8524/6,7/calexp-HSC-G-8524-6,7.fits']

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
