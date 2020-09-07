#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 14:38:27 2020

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

def find_loc_max(data, neighborhood_size = 8, threshold = 5):
    neighborhood_size = neighborhood_size
    threshold = threshold
    data_max = filters.maximum_filter(data, neighborhood_size) 
    maxima = (data == data_max)
    data_min = filters.minimum_filter(data, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0
    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
    x, y = [], []
    for dy,dx in slices:
        x_center = (dx.start + dx.stop - 1)/2
        x.append(x_center)
        y_center = (dy.start + dy.stop - 1)/2    
        y.append(y_center)
    return x, y

def search_local_max(img, radius=100, view=False):
    from .cutout_tools import cutout
    from .astro_tools import plt_fits
    PSFx, PSFy =find_loc_max(img)
    PSF_locs = []
    ct = 0
    for i in range(len(PSFx)):
        cut_img = cutout(img, [PSFx[i], PSFy[i]], radius=radius)
        cut_img[np.isnan(cut_img)] = 0
        if np.sum(cut_img==0)  > len(cut_img)**2/5:
            continue
        PSF_locs.append([PSFx[i], PSFy[i]])
        if view ==True:
            print("plot for position: [{0}, {1}]".format(PSFx[i], PSFy[i]), "idx:", ct)
            print("total flux:", cut_img.sum())
            print("measure FWHM:", np.round(measure_FWHM(cut_img),3))
            plt_fits(cut_img)
            print("================")
            ct += 1
    return  PSF_locs
    

def measure_FWHM(image, measure_radius = 10):
    seed_num = 2*measure_radius+1
    frm = len(image)
    q_frm = int(frm/4)
    x_center = np.where(image == image[q_frm:-q_frm,q_frm:-q_frm].max())[1][0]
    y_center = np.where(image == image[q_frm:-q_frm,q_frm:-q_frm].max())[0][0]
    
    x_n = np.asarray([image[y_center][x_center+i] for i in range(-measure_radius, measure_radius+1)]) # The x value, vertcial 
    y_n = np.asarray([image[y_center+i][x_center] for i in range(-measure_radius, measure_radius+1)]) # The y value, horizontal 
    xy_n = np.asarray([image[y_center+i][x_center+i] for i in range(-measure_radius, measure_radius+1)]) # The up right value, horizontal
    xy__n =  np.asarray([image[y_center-i][x_center+i] for i in range(-measure_radius, measure_radius+1)]) # The up right value, horizontal
    from astropy.modeling import models, fitting
    g_init = models.Gaussian1D(amplitude=y_n.max(), mean=measure_radius, stddev=1.5)
    fit_g = fitting.LevMarLSQFitter()
    g_x = fit_g(g_init, range(seed_num), x_n)
    g_y = fit_g(g_init, range(seed_num), y_n)
    g_xy = fit_g(g_init, range(seed_num), xy_n)
    g_xy_ = fit_g(g_init, range(seed_num), xy__n)
    
    FWHM_ver = g_x.stddev.value * 2.355  # The FWHM = 2*np.sqrt(2*np.log(2)) * stdd = 2.355*stdd
    FWHM_hor = g_y.stddev.value * 2.355
    FWHM_xy = g_xy.stddev.value * 2.355 * np.sqrt(2.)
    FWHM_xy_ = g_xy_.stddev.value * 2.355 * np.sqrt(2.)
    return FWHM_ver, FWHM_hor, FWHM_xy, FWHM_xy_
