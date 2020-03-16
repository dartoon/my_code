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
            
name='NGC1376'
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
data_fix=changebox(data_fix,center=[242,255], half_size=[6, 9], shift_center=[258, 243])
data_fix=changebox(data_fix,center=[42,252], half_size=[6, 7], shift_center=[29, 235])
data_fix=changebox(data_fix,center=[294,129], half_size=[5, 8], shift_center=[293, 105])
data_fix=changebox(data_fix,center=[281,204], half_size=[1, 1], shift_center=[283, 200])
data_fix=changebox(data_fix,center=[35,275], half_size=[1, 1], shift_center=[31, 272])
data_fix=changebox(data_fix,center=[18,279], half_size=[1, 1], shift_center=[17, 275])
#==============================================================================
# if need smooth
#==============================================================================
data_fix=smooth_edge_bump(data_fix, r_out=65, thre=0.5)
data_fix=smooth_edge_bump(data_fix, r_out=128, thre=0.08)


data_fix -= 0.038 # np.median(data_fix[0:35,0:35])  # subtract a "background" value; or -np.median(data_rehalf_sized[0:10,0:10])

plt.matshow(np.log10(data_fix.T),origin='lower')
plt.colorbar()
plt.show()

pyfits.PrimaryHDU(data_fix.T).writeto('fix/{0}_fix.fits'.format(name),clobber=True)

from shutil import copyfile
copyfile('cosmetology.py', 'fix/cosmetology_{0}.py'.format(name))
