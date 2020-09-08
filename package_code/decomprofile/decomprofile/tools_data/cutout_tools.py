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
        radius: the cutout box size.
    Return
    --------
        A cutout image stamp, frame size in odd number.
    """
    region = pix_region(center, radius=radius)
    cut = region.to_mask(mode='exact')
    cut_image = cut.cutout(image)
    return cut_image

def cut_center_auto(image, center, radius, kernel = 'center_bright', return_center=False,
                      if_plot=False):
    """
    Automaticlly cutout out a image, so that the central pixel is either the "center_bright" or "center_gaussian".
    
    Parameter
    --------
        image: large frame 2D image data;
        center: the center position to cutout;
        kernel: the way to define the central pixel, either 
            'center_bright': 
                Cutout at the brightest pixel as center
            or 
            'center_gaussian': 
                Cutout at the Gaussian center
        radius: the cutout box size;
        return_center: if return the finally used center value;
        if_plot: if plot the zoom in center of the cutout stamp.
    
    Warning  
    --------
    Frame size shouldn't be too larger that exceed the target's frame otherwise the Gaussion fitting and the max pixel could miss targets 
        
    Return
    --------
        A cutout image; (if return_center= True, a central pixel position used for the cutout would be returned)
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
        center_pos = (temp_center.astype(int) + np.round(center_shift))
        cutout_image = cutout(image=image, center=center_pos, radius=radius)
        plt_center = img_init[frm_q:-frm_q,frm_q:-frm_q].shape
        if if_plot==True:
            plt.plot(plt_center[0]/2-0.5, plt_center[1]/2-0.5, color='c', marker='+', ms=ms, mew=mew)  #-0.5 to shift the "+" to the center
            plt.imshow(cutout_image[frm_q:-frm_q,frm_q:-frm_q], origin='lower', norm=LogNorm())
            plt.show()
    elif kernel == 'center_gaussian':
        # test_center = frm_q + centroid_2dg(img_init[frm_q:-frm_q,frm_q:-frm_q])
        gauss_center = centroid_2dg(img_init)
        center_shift = gauss_center - radius
        center = np.asarray(center)
        center_pos = center.astype(int) + center_shift
        img_init = cutout(image=image, center=center_pos, radius=radius)
        if if_plot==True :
            fig, ax = plt.subplots(1, 1)
            plt_center = img_init[frm_q:-frm_q,frm_q:-frm_q].shape
            plt.plot(plt_center[0]/2-0.5, plt_center[1]/2-0.5, color='r', marker='+', ms=ms, mew=mew)
            plt.imshow(img_init[frm_q:-frm_q,frm_q:-frm_q], origin='lower', norm=LogNorm())
            plt.show()
        cutout_image = img_init
    else:
        raise ValueError("kernel is not defined")
    if return_center==False:
        return cutout_image
    elif return_center==True:
        return cutout_image, center_pos

def plot_overview(img, center_QSO, c_psf_list=None, label=None, ifsave=False):
    """
    Plot the overview of the image, highlight the location of the QSO and PSFs
    
    Parameter
    --------
        img: A FOV image
        center_QSO: The central position of the pixels of the QSO
        c_psf_list: A list of PSF positions.
        label: define label if want to lable this plot
        
        
    Return
    --------
        No returns
    """    
    import copy, matplotlib
    my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
    my_cmap.set_bad('black')
    vmax = 2.2
    vmin = 1.e-2
    QSO_box_size = np.min(img.shape)/72
    PSF_box_size = np.min(img.shape)/109
    fig = plt.figure(figsize=(15,15))
    ax=fig.add_subplot(1,1,1)
    ax.imshow(img,origin='lower', cmap=my_cmap, norm=LogNorm(), vmin=vmin, vmax=vmax)
    QSO_reg = pix_region(center_QSO, radius= QSO_box_size)
    QSO_mask = QSO_reg.to_mask(mode='center')
    ax.text(center_QSO[0]-2*QSO_box_size, center_QSO[1]+1.5*QSO_box_size, 'QSO',color='white', fontsize=20)
    ax.add_patch(QSO_mask.bbox.as_artist(facecolor='none', edgecolor='white', linewidth=2))
    name = 'PSF'
    count=0
    if c_psf_list is not None:
        for i in range(len(c_psf_list)):
            PSF_reg = pix_region(c_psf_list[i], radius= PSF_box_size)
            PSF_mask = PSF_reg.to_mask(mode='center')
            ax.add_patch(PSF_mask.bbox.as_artist(facecolor='none', edgecolor='blue', linewidth=2))
            ax.text(c_psf_list[i][0]-2*PSF_box_size, c_psf_list[i][1]+2*PSF_box_size, '{1}{0}'.format(count, name),color='white', fontsize=15)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            count += 1
#    plt.colorbar(cax)
    if not label == None:
        ax.text(len(img)*0.05, len(img)*0.8, label,color='white', fontsize=30)
    if ifsave == True:
        fig.savefig('QSO_{0}_loc.pdf'.format(name))    

def string_find_between(s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

def cr_mask(image, filename='test_circle.reg'):
    '''
    The creat a mask using a DS9 .reg file. The pixels in the region are 0, ouside ones are 1.
    
    Parameter
    --------
        image: a 2D array image as a template frame.
        filename: full name of the region file.
        
    Return
    --------
        A image.shape array. Pixels in the region is 0, otherwise 1.
    '''
    with open(filename, 'r') as input_file:
        reg_string=input_file.read().replace('\n', '')
    if "physicalcircle" in reg_string:
        abc=string_find_between(reg_string, "(", ")")
        reg_info=np.fromstring(abc, dtype=float, sep=',')
        center, radius = reg_info[:2]-1 , reg_info[2]
        region = pix_region(center, radius)
        box = 1-region.to_mask(mode='center').data
    elif "physicalbox" in reg_string:
        abc=string_find_between(reg_string, "(", ")")
        reg_info=np.fromstring(abc, dtype=float, sep=',')
        center = reg_info[:2] - 1
        x_r, y_r = reg_info[2:4]  # x_r is the length of the x, y_r is the length of the y
        box = np.zeros([np.int(x_r)+1, np.int(y_r)+1]).T
    else:
        print(reg_string)
        raise ValueError("The input reg is un-defined yet")
    frame_size = image.shape
    box_size = box.shape
    x_edge = np.int(center[1]-box_size[0]/2) #The position of the center is x-y switched.
    y_edge = np.int(center[0]-box_size[1]/2)
    mask = np.ones(frame_size)
    mask_box_part = mask[x_edge:x_edge+box_size[0],y_edge: y_edge + box_size[1]]
    mask_box_part *= box
    return mask