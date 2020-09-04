#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 15:12:15 2020

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from regions import PixCoord, CirclePixelRegion 
from matplotlib.colors import LogNorm

def pix_region(center=[49.0,49.0], radius=5):
    """
    Creat a region file, in pixel units
    
    Parameter
    --------
        center: The center of the region, with [reg_x, reg_y];
        radius: The radius of the region.
    Return
    --------
        A region which is ds9-like.
    """
    center= PixCoord(x=center[0],y=center[1])
    region = CirclePixelRegion(center, radius)
    return region

def cutout(image, center, radius):
    """
    Cutout a stamp image from a large frame image
    
    Parameter
    --------
        image: large frame 2D image data;
        center: the center position to cutout;
        radius: the cutout box size
    Return
    --------
        A cutout image stamp
    """
    region = pix_region(center, radius=radius)
    cut = region.to_mask(mode='exact')
    cut_image = cut.cutout(image)
    return cut_image

def cut_center_bright(image, center, radius, kernel = 'center_bright', return_center=False,
                      if_plot=True):
    """
    Automaticlly cutout out a image, so that the central pixel is either the "center_bright" or "center_gaussian".
    
    Parameter
    --------
        image: large frame 2D image data;
        center: the center position to cutout;
        kernel: the way to define the central pixel, either 
            'center_bright'
            or 
            'center_gaussian'
        radius: the cutout box size
        return_center: if return the finally used center value
        if_plot: if plot the cutout
        
    Return
    --------
        A sth sth
    """
    from photutils import centroid_2dg
    temp_center = np.asarray(center)
#    print temp_center.astype(int)
    radius = radius
    img_init = cutout(image=image, center=temp_center.astype(int), radius=radius)
    frm_q = int(len(img_init)/2.5)  #Aa quarter scale of the frame
    ms, mew = 30, 2.
    if kernel == 'center_bright':
        test_center =  np.asarray(np.where(img_init == img_init[frm_q:-frm_q,frm_q:-frm_q].max()))[:,0]
        center_shift = np.array((test_center- radius))[::-1]
        center = (temp_center.astype(int) + np.round(center_shift))
        cutout_image = cutout(image=image, center=center, radius=radius)
        plt_center = img_init[frm_q:-frm_q,frm_q:-frm_q].shape
        if if_plot==True:
            marker = '+'
            plt.plot(plt_center[0]/2, plt_center[1]/2, color='c', marker='+', ms=ms, mew=mew)
            plt.imshow(cutout_image[frm_q:-frm_q,frm_q:-frm_q], origin='lower')
            plt.show()
    elif kernel == 'center_gaussian':
        for i in range(3):
            test_center = frm_q + centroid_2dg(img_init[frm_q:-frm_q,frm_q:-frm_q])
            if i ==2 and if_plot==True :
                print(test_center)
                fig, ax = plt.subplots(1, 1)
                ax.imshow(img_init[frm_q:-frm_q,frm_q:-frm_q], origin='lower', norm=LogNorm())
                marker = '+'
                plt.plot(test_center[0]-frm_q, test_center[1]-frm_q, color='b', marker=marker, ms=ms, mew=mew)
                plt.show()
                print('test_center - radius', test_center, radius)
            center_shift = np.array((test_center - radius))
            if i ==2 and if_plot==True :
                print('center_shift',center_shift)
            center = np.asarray(center)
            center = center.astype(int) + center_shift
            img_init = cutout(image=image, center=center, radius=radius)
            if i ==2 and if_plot==True :
                plt_center = img_init[frm_q:-frm_q,frm_q:-frm_q].shape
                plt.plot(plt_center[0]/2, plt_center[1]/2, color='r', marker=marker, ms=ms, mew=mew)
                plt.imshow(img_init[frm_q:-frm_q,frm_q:-frm_q], origin='lower', norm=LogNorm())
                plt.show()
        cutout_image = img_init
    else:
        raise ValueError("kernel is not defined")
    if return_center==False:
        return cutout_image
    elif return_center==True:
        return cutout_image, center
