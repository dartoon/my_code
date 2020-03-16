#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:15:07 2017

@author: dxh

data_fix -= np.median(data_fix[0:35,0:35])  # subtract a "background" value; or -np.median(data_resized[0:10,0:10])
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
name='NGC7742'

def changebox(data, center, size, shift_center):
    '''
    data=np.array like
    center= [x_center, y_center]
    size= [x_range, y_rang]  --- total size is 2*x_range+1
    shift = [x_shift, y_shift]
    '''
    start_point=[center[0]-size[0], center[1]-size[1]]
    for i in range(size[0]*2+1):
        for j in range(size[1]*2+1):
            data[start_point[0]+i,start_point[1]+j]=data[shift_center[0]+i-size[0],shift_center[1]+j-size[1]]
    return data
    

hdu = pyfits.open('../resized_source/{0}_resized.fits'.format(name))
data = hdu[0].data.T
hdu.close()

plt.matshow(np.log10(data.T),origin='lower')
plt.colorbar()
plt.show()
data_fix=data
data_fix=changebox(data_fix,center=[231,359], size=[14, 13], shift_center=[276, 340])
data_fix=changebox(data_fix,center=[251,37], size=[11, 10], shift_center=[220, 30])

plt.matshow(np.log10(data_fix.T),origin='lower')
plt.colorbar()
plt.show()

data_fix -= 0.0823062  # subtract a "background" value; or -np.median(data_resized[0:10,0:10])
pyfits.PrimaryHDU(data_fix.T).writeto('{0}_fix.fits'.format(name),clobber=True)

