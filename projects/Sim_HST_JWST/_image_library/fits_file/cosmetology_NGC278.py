#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:15:07 2017

@author: dxh

data_fix -= np.median(data_fix[0:35,0:35])  # subtract a "background" value; or -np.median(data_rehalf_sized[0:10,0:10])
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import sys
sys.path.insert(0,'../../evil_code/')
from changebox import changebox
from smooth_edge_bump import smooth_edge_bump
            
name='NGC278'
hdu = pyfits.open('{0}_resized.fits'.format(name))
data = hdu[0].data.T
hdu.close()

plt.matshow(np.log10(data.T),origin='lower')
plt.colorbar()
plt.show()

data_fix=data
data_fix=changebox(data_fix,center=[102,228], half_size=[5/2, 5/2], shift_center=[117, 235])
data_fix=changebox(data_fix,center=[233,201], half_size=[9/2, 9/2], shift_center=[226, 214])
#
data_fix=smooth_edge_bump(data_fix, r_out=85, thre=0.9)
data_fix=smooth_edge_bump(data_fix, r_out=115, thre=0.4)

plt.matshow(np.log10(data_fix.T),origin='lower')
plt.colorbar()
plt.show()

data_fix -= 0.0620392 # np.median(data_fix[0:35,0:35])  # subtract a "background" value; or -np.median(data_rehalf_sized[0:10,0:10])
pyfits.PrimaryHDU(data_fix.T).writeto('fix/{0}_fix.fits'.format(name),clobber=True)

from shutil import copyfile
copyfile('cosmetology.py', 'fix/cosmetology_{0}.py'.format(name))
