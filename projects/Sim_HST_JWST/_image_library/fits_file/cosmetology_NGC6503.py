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
            
name='NGC6503'
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
data_fix=changebox(data_fix,center=[65,250], half_size=[7, 6], shift_center=[83, 259])
data_fix=changebox(data_fix,center=[106,204], half_size=[7, 6], shift_center=[118, 209])
data_fix=changebox(data_fix,center=[235,225], half_size=[3, 3], shift_center=[244, 225])
data_fix=changebox(data_fix,center=[181,256], half_size=[3, 3], shift_center=[189, 258])
data_fix=changebox(data_fix,center=[14,95], half_size=[4, 5], shift_center=[20, 111])
data_fix=changebox(data_fix,center=[242,88], half_size=[3, 3], shift_center=[234, 79])
data_fix=changebox(data_fix,center=[10,196], half_size=[2, 3], shift_center=[18, 201])
#==============================================================================
# if need smooth
#==============================================================================
data_fix=smooth_edge_bump(data_fix, r_out=95, thre=1.1)


#data_fix -= 0.0965518# np.median(data_fix[0:35,0:35])  # subtract a "background" value; or -np.median(data_rehalf_sized[0:10,0:10])

plt.matshow(np.log10(data_fix.T),origin='lower')
plt.colorbar()
plt.show()

pyfits.PrimaryHDU(data_fix.T).writeto('fix/{0}_fix.fits'.format(name),clobber=True)

from shutil import copyfile
copyfile('cosmetology.py', 'fix/cosmetology_{0}.py'.format(name))
