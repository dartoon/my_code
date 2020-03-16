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
            
name='NGC5584'
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
data_fix=changebox(data_fix,center=[282,158], half_size=[6, 9], shift_center=[281, 137])
data_fix=changebox(data_fix,center=[235,8], half_size=[4, 4], shift_center=[238, 18])
data_fix=changebox(data_fix,center=[170,278], half_size=[6, 5], shift_center=[154, 277])
data_fix=changebox(data_fix,center=[95,187], half_size=[3, 4], shift_center=[103, 186])
data_fix=changebox(data_fix,center=[159,241], half_size=[4, 4], shift_center=[171, 254])
data_fix=changebox(data_fix,center=[200,250], half_size=[4, 4], shift_center=[183, 244])
data_fix=changebox(data_fix,center=[295,177], half_size=[4, 5], shift_center=[295, 191])
data_fix=changebox(data_fix,center=[295,211], half_size=[4, 5], shift_center=[295, 195])
data_fix=changebox(data_fix,center=[44,255], half_size=[3, 4], shift_center=[50, 264])

#==============================================================================
# if need smooth
#==============================================================================
data_fix=smooth_edge_bump(data_fix, r_out=17, thre=0.22)
data_fix=smooth_edge_bump(data_fix, r_out=115, thre=0.18)
#data_fix -= 0.0965518 # np.median(data_fix[0:35,0:35])  # subtract a "background" value; or -np.median(data_rehalf_sized[0:10,0:10])



data_fix[data_fix!=0] -= 0.035  # subtract a "background" value; or -np.median(data_resized[0:10,0:10])
data_fix[data_fix<-0.008] = 0 # subtract a "background" value; or -np.median(data_resized[0:10,0:10])

plt.matshow(np.log10(data_fix.T),origin='lower')
plt.colorbar()
plt.show()

pyfits.PrimaryHDU(data_fix.T).writeto('fix/{0}_fix.fits'.format(name),clobber=True)

from shutil import copyfile
copyfile('cosmetology.py', 'fix/cosmetology_{0}.py'.format(name))
