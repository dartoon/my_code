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
            
name='NGC2397'
hdu = pyfits.open('{0}_resized.fits'.format(name))
data = hdu[0].data.T
hdu.close()

plt.matshow(np.log10(data.T),origin='lower')
plt.colorbar()
plt.show()

data_fix=data
data_fix=changebox(data_fix,center=[95,346], half_size=[7/2, 7/2], shift_center=[75, 336])
data_fix=changebox(data_fix,center=[112,272], half_size=[9/2, 9/2], shift_center=[119, 279])
data_fix=changebox(data_fix,center=[159,16], half_size=[10/2, 8/2], shift_center=[181,22])

data_fix=smooth_edge_bump(data_fix, r_out=27, thre=2.45)
data_fix=smooth_edge_bump(data_fix, r_out=59, thre=1.8)
data_fix=smooth_edge_bump(data_fix, r_out=94, thre=1.12)
data_fix=smooth_edge_bump(data_fix, r_out=118, thre=0.3)
data_fix=smooth_edge_bump(data_fix, r_out=126, thre=0.2)

plt.matshow(np.log10(data_fix.T),origin='lower')
plt.colorbar()
plt.show()

data_fix -= 0.00558646 # np.median(data_fix[0:35,0:35])  # subtract a "background" value; or -np.median(data_rehalf_sized[0:10,0:10])
pyfits.PrimaryHDU(data_fix.T).writeto('fix/{0}_fix.fits'.format(name),clobber=True)

from shutil import copyfile
copyfile('cosmetology.py', 'fix/cosmetology_{0}.py'.format(name))
