#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 14:39:05 2018

@author: dartoon

Test create mask auto
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

from mask_objects import mask_obj

import os
path = os.getcwd()
ID = path.split('/')[-2]

QSO_img  = pyfits.getdata('test.fits')
masks  = mask_obj(img=QSO_img, exp_sz=1.2)

obj_mask = masks[0] + masks[1] + masks[2]
obj_mask = [obj_mask==0][0] + 0

plt.imshow(obj_mask)
print "the yellow area(==1) is the area to fit"

                      
#kwargs_numerics['mask'] = obj_mask