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
            
name='NGC1309'
hdu = pyfits.open('resized_source/{0}_resized.fits'.format(name))
data = hdu[0].data.T
hdu.close()

plt.matshow(np.log10(data.T),origin='lower')
plt.colorbar()
plt.show()

data_fix=data
#==============================================================================
# if changebox
#==============================================================================
data_fix=changebox(data_fix,center=[153,235], half_size=[10, 9], shift_center=[156, 76],y_re=-1)
data_fix=changebox(data_fix,center=[285,86], half_size=[11, 10], shift_center=[272, 56])
data_fix=changebox(data_fix,center=[283,277], half_size=[6, 6], shift_center=[259, 281])
data_fix=changebox(data_fix,center=[291,33], half_size=[8, 9], shift_center=[257, 18])
#==============================================================================
# if need smooth
#==============================================================================
data_fix=smooth_edge_bump(data_fix, r_out=25, thre=2.6)

data_fix -= 0.0032778-0.001 # np.median(data_fix[0:35,0:35])  # subtract a "background" value; or -np.median(data_rehalf_sized[0:10,0:10])

plt.matshow(np.log10(data_fix.T),origin='lower')
plt.colorbar()
plt.show()

pyfits.PrimaryHDU(data_fix.T).writeto('fix/{0}_fix.fits'.format(name),clobber=True)

from shutil import copyfile
copyfile('cosmetology.py', 'fix/cosmetology_{0}.py'.format(name))
