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
            
name='ESO498G5'
hdu = pyfits.open('{0}_resized.fits'.format(name))
data = hdu[0].data.T
hdu.close()

plt.matshow(np.log10(data.T),origin='lower')
plt.colorbar()
plt.show()

data_fix=data
#==============================================================================
# if changebox
#==============================================================================
data_fix=changebox(data_fix,center=[99,269], half_size=[15/2, 17/2], shift_center=[119, 291])
data_fix=changebox(data_fix,center=[255,298], half_size=[22/2, 21/2], shift_center=[276, 263])
data_fix=changebox(data_fix,center=[197,231], half_size=[7/2, 6/2], shift_center=[204, 225])
data_fix=changebox(data_fix,center=[215,219], half_size=[7/2, 6/2], shift_center=[218, 211])
data_fix=changebox(data_fix,center=[243,75], half_size=[8/2, 7/2], shift_center=[229, 65])
data_fix=changebox(data_fix,center=[291,51], half_size=[8/2, 8/2], shift_center=[284, 41])
data_fix=changebox(data_fix,center=[87,61], half_size=[7/2, 7/2], shift_center=[98, 60])
data_fix=changebox(data_fix,center=[232,18], half_size=[6/2, 6/2], shift_center=[243, 26])
data_fix=changebox(data_fix,center=[322,72], half_size=[10/2, 18/2], shift_center=[311, 58])
data_fix=changebox(data_fix,center=[309,189], half_size=[20/2, 16/2], shift_center=[308, 164])
#==============================================================================
# if need smooth
#==============================================================================
data_fix=smooth_edge_bump(data_fix, r_out=27, thre=0.4)

data_fix=smooth_edge_bump(data_fix, r_out=114, thre=0.135)
data_fix=smooth_edge_bump(data_fix, r_out=162, thre=0.11)

plt.matshow(np.log10(data_fix.T),origin='lower')
plt.colorbar()
plt.show()

data_fix -= 0.0115312 # np.median(data_fix[0:35,0:35])  # subtract a "background" value; or -np.median(data_rehalf_sized[0:10,0:10])
pyfits.PrimaryHDU(data_fix.T).writeto('fix/{0}_fix.fits'.format(name),clobber=True)

from shutil import copyfile
copyfile('cosmetology.py', 'fix/cosmetology_{0}.py'.format(name))
