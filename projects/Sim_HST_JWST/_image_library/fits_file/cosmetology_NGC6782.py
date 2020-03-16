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
            
name='NGC6782'
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
data_fix=changebox(data_fix,center=[76,238], half_size=[18, 18], shift_center=[101, 279])
data_fix=changebox(data_fix,center=[47,276], half_size=[18, 18], shift_center=[78, 313])
data_fix=changebox(data_fix,center=[260,356], half_size=[10, 9], shift_center=[282, 328])
data_fix=changebox(data_fix,center=[354,249], half_size=[9, 9], shift_center=[363, 225])
data_fix=changebox(data_fix,center=[178,9], half_size=[6, 4], shift_center=[201, 9])
data_fix=changebox(data_fix,center=[292,338], half_size=[6, 4], shift_center=[308, 318])

#==============================================================================
# if need smooth
#==============================================================================
data_fix=smooth_edge_bump(data_fix, r_out=41.6352, thre=1.76)
data_fix=smooth_edge_bump(data_fix, r_out=137.354, thre=0.89)

data_fix -= 0.18# np.median(data_fix[0:35,0:35])  # subtract a "background" value; or -np.median(data_rehalf_sized[0:10,0:10])

plt.matshow(np.log10(data_fix.T),origin='lower')
plt.colorbar()
plt.show()

pyfits.PrimaryHDU(data_fix.T).writeto('fix/{0}_fix.fits'.format(name),clobber=True)

from shutil import copyfile
copyfile('cosmetology.py', 'fix/cosmetology_{0}.py'.format(name))
