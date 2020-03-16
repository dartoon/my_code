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
            
name='NGC1084'
hdu = pyfits.open('{0}_resized.fits'.format(name))
data = hdu[0].data.T
hdu.close()

plt.matshow(np.log10(data.T),origin='lower')
plt.colorbar()
plt.show()

data_fix=data
data_fix=changebox(data_fix,center=[72,205], half_size=[6/2, 6/2], shift_center=[79, 198])
data_fix=changebox(data_fix,center=[81,167], half_size=[6/2, 6/2], shift_center=[59, 158])
data_fix=changebox(data_fix,center=[335,316], half_size=[8/2, 7/2], shift_center=[345, 306])
data_fix=changebox(data_fix,center=[380,255], half_size=[6/2, 8/2], shift_center=[382, 242])
data_fix=changebox(data_fix,center=[99,112], half_size=[9/2, 9/2], shift_center=[111, 103])
data_fix=changebox(data_fix,center=[330,172], half_size=[5/2, 5/2], shift_center=[336, 165])

data_fix=changebox(data_fix,center=[245,370], half_size=[11/2, 11/2], shift_center=[261, 359])
data_fix=changebox(data_fix,center=[244,33], half_size=[33/2, 25/2], shift_center=[215, 35])

#data_fix=smooth_edge_bump(data_fix, r_out=123, thre=0.9)
#data_fix=smooth_edge_bump(data_fix, r_out=115, thre=0.4)

plt.matshow(np.log10(data_fix.T),origin='lower')
plt.colorbar()
plt.show()

data_fix -= 0.00876185 # np.median(data_fix[0:35,0:35])  # subtract a "background" value; or -np.median(data_rehalf_sized[0:10,0:10])
pyfits.PrimaryHDU(data_fix.T).writeto('fix/{0}_fix.fits'.format(name),clobber=True)

from shutil import copyfile
copyfile('cosmetology.py', 'fix/cosmetology_{0}.py'.format(name))
