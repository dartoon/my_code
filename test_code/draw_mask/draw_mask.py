#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 22:43:59 2018

@author: Dartoon

plot mask for one image
"""

from regions import PixCoord, CirclePixelRegion 
import numpy as np
import matplotlib.pylab as plt
import astropy.io.fits as pyfits
import glob

def pix_region(center=([49,49]), radius=5):
    '''
    Creat a region file, with pixel units
    
    Parameter
    --------
        center: The center of the region, with ([reg_x, reg_y]);
        radius: The radius of the region.
        
    Return
    --------
        A region which is ds9-like.
    '''
    center= PixCoord(x=center[0],y=center[1])
    region = CirclePixelRegion(center, radius)
    return region

def cr_mask(image, filename='test_circle.reg'):
    '''
    The creat a mask with a .reg file. The pixels in the region is 0, otherwise 1.
    
    Parameter
    --------
        filename: filename of the .reg
        image:
        
    Return
    --------
        A image.shape array. Pixels in the region is 0, otherwise 1.
    '''
    ##### Note the numpy starts from 0, especially in the center,
    ####!!!Need to check the center, shape of the boxtype, the size?
    with open(filename, 'r') as input_file:
        reg_string=input_file.read().replace('\n', '')
    if "physicalcircle" in reg_string:
        abc=string_find_between(reg_string, "(", ")")
        reg_info=np.fromstring(abc, dtype=float, sep=',')
        center, radius = reg_info[:2]-1, reg_info[2]
        region = pix_region(center, radius)
        box = 1-region.to_mask(mode='center').data
    elif "physicalbox" in reg_string:
        abc=string_find_between(reg_string, "(", ")")
        reg_info=np.fromstring(abc, dtype=float, sep=',')
        center = reg_info[:2] - 1
        x_r, y_r = reg_info[2:4]  # x_r is the length of the x, y_r is the length of the y
        box = np.zeros([np.int(x_r)+1, np.int(y_r)+1]).T
    else:
        print reg_string
        raise ValueError("The input reg is un-defined yet")
    frame_size = image.shape
    box_size = box.shape
    x_edge = np.int(center[1]-box_size[0]/2) #The position of the center is x-y switched.
    y_edge = np.int(center[0]-box_size[1]/2)
    mask = np.ones(frame_size)
    mask_box_part = mask[x_edge:x_edge+box_size[0],y_edge: y_edge + box_size[1]]
    mask_box_part *= box
    return mask

def string_find_between(s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


mask_list = glob.glob('mask_reg*.reg')
data = pyfits.getdata('test_img.fits')

mask=np.ones_like(data)
          
for i in range(len(mask_list)):
    mask *= cr_mask(image=data, filename=mask_list[i])

mask = 1-mask
plt.imshow(mask, origin='low')
plt.show()
pyfits.PrimaryHDU(mask).writeto('mask.fits',overwrite=True)