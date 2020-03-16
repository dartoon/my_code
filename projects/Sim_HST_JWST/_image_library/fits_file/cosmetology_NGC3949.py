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
name='NGC3949'

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
    

hdu = pyfits.open('{0}_resized.fits'.format(name))
data = hdu[0].data.T
hdu.close()
data -= np.median(data[20:,0:20])

plt.matshow(np.log10(data.T),origin='lower')
plt.colorbar()
plt.show()
data_fix=data
data_fix=changebox(data_fix,center=[304,160], size=[3, 3], shift_center=[296, 148])
data_fix=changebox(data_fix,center=[118,396], size=[8, 3], shift_center=[135, 396])
data_fix=changebox(data_fix,center=[134,10], size=[3, 4], shift_center=[146, 13])

p1_l=[307.,399.]
p2_l=[346.,344.]
a1=(p2_l[1]-p1_l[1])/(p2_l[0]-p1_l[0])
b1=p1_l[1]-p1_l[0]*a1

p1_r=[346.,344.]
p2_r=[400.,375.]
a2=(p2_r[1]-p1_r[1])/(p2_r[0]-p1_r[0])
b2=p1_r[1]-p1_r[0]*a2
for i in range(306,len(data)):
    for j in range(353,len(data)):
        if (j-a1*i-b1)>0 and (j-a2*i-b2)>0:
            data_fix[i,j]=data_fix[len(data)-i,len(data)-j]
#            print i, j 
blank=np.where((data_fix<-0.01))
for i in range(len(blank[0])):
      data_fix[blank[0][i],blank[1][i]]=data_fix[(len(data)-blank[0][i]-1),(len(data)-blank[1][i]-1)]

data_fix -= -0.038  # subtract a "background" value; or -np.median(data_resized[0:10,0:10])
plt.matshow(np.log10(data_fix.T),origin='lower')
plt.colorbar()
plt.show()

pyfits.PrimaryHDU(data_fix.T).writeto('fix/{0}_fix.fits'.format(name),clobber=True)

from shutil import copyfile
copyfile('cosmetology.py', 'fix/cosmetology_{0}.py'.format(name))
