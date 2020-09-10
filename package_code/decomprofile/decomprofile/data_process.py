#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 09:10:39 2020

@author: Xuheng Ding

A class to process the data
"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from astropy.wcs import WCS
from decomprofile.tools_data.measure_tools import measure_bkg
from decomprofile.tools_data.cutout_tools import cut_center_auto, cutout
from copy import deepcopy
from matplotlib.colors import LogNorm
from decomprofile.tools_data.astro_tools import plt_fits, read_pixel_scale

class DataProcess(object):
    """
    A class to Process the data, including the following feature:
        - automaticlly estimate and remove background light.
        - cutout the target photo stamp.
        - search all the avaiable PSF in the field.
        - creat mask for the objects.
        - measure the target surface brightness profile, PSF FWHM, background.
    """
    def __init__(self, fov_image, target_pos, pos_type = 'pixel', header=None, exptime = None,
                 rm_bkglight = True, if_plot = False, nsigma=2, npixels=25, dilate_size=11):
        """
        Parameter
        --------
            data_image: 2D array
            The field of view image of the data.
            
            target_pos: list or tuple or array, length = 2
            The position of the target.
            
            pos_type: string, 'pixel' or 'wcs'
            Define the position of the target, i.e., if the position is in 'pixel' or 'wcs'.
                
            header: io.fits.header
            The header information given by the fits file. 
            Note: should including the exposure time and WCS information.
            
            exptime: float or 2D array
            The exposure time of the data in (s) a the exptime_map
            
        """
        if pos_type == 'pixel':
            self.target_pos = target_pos
        elif pos_type == 'wcs':
            wcs = WCS(header)
            self.target_pos = wcs.all_world2pix([[target_pos[0], target_pos[1]]], 1)
        else:
            raise ValueError("'pos_type' is should be either 'pixel' or 'wcs'.")

        self.exptime = exptime
        self.if_plot = if_plot    
        self.header = header
        self.deltaPix = read_pixel_scale(header)
        if rm_bkglight == True:
            bkglight = measure_bkg(fov_image, if_plot=if_plot, nsigma=nsigma,
                                   npixels=npixels,dilate_size=dilate_size)
            fov_image = fov_image-bkglight
        self.fov_image = fov_image
        


    def generate_target_materials(self, cut_kernel = 'center_gaussian',  target_radius=60, 
                                  bkg_std = None, create_mask = False, if_plot=None, **kwargs):
        """
        Produce the materials that would be used for the fitting.
        
        Parameter
        --------
            target_radius: int or float
            The radius of aperture to cutout the target
            
            bkg_std: The blash of blash
            
        """
        if if_plot == None:
            if_plot = self.if_plot
        
        target_stamp, target_cut_pos = cut_center_auto(image=self.fov_image, center= self.target_pos, 
                                          kernel = cut_kernel, radius=target_radius,
                                          return_center=True, if_plot=if_plot)
        if bkg_std == None:
            from decomprofile.tools_data.measure_tools import esti_bgkstd
            target_2xlarger_stamp = cutout(image=self.fov_image, center= target_cut_pos, radius=target_radius*2)
            self.bkg_std = esti_bgkstd(target_2xlarger_stamp, if_plot=if_plot)
        exptime = deepcopy(self.exptime)
        if exptime is None:
            if 'EXPTIME' in self.header.keys():
                exptime = self.header['EXPTIME']
            else:
                raise ValueError("No Exposure time information in the header, should input a value.")
        if isinstance(exptime, np.ndarray):
            exptime_stamp = cutout(image=self.exptime, center= target_cut_pos, radius=target_radius)
        noise_map = np.sqrt(abs(target_stamp/exptime_stamp) + self.bkg_std**2)
        target_mask = np.ones_like(target_stamp)
        if create_mask == True:
            from decomprofile.tools_data.measure_tools import detect_obj, mask_obj
            apertures = detect_obj(target_stamp, if_plot=True, **kwargs)
            select_idx = input('Input directly the a obj idx to mask, use space between each id:\n')
            select_idx = select_idx.split(" ")
            select_idx = [int(select_idx[i]) for i in range(len(select_idx)) if select_idx[i].isnumeric()]
            apertures = [apertures[i] for i in select_idx]
            mask_list = mask_obj(target_stamp, apertures, if_plot=False)
            for i in range(len(mask_list)):
                target_mask *= mask_list[i]
        if if_plot:
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 10))
            ax1.imshow(target_stamp, origin='lower', norm=LogNorm())
            ax1.set_title('Cutout target')
            ax2.imshow(noise_map, origin='lower', norm=LogNorm())
            ax2.set_title('Noise map')
            ax3.imshow(target_stamp * target_mask, origin='lower', norm=LogNorm())
            ax3.set_title('data * mask')
            plt.show()  
        self.target_stamp = target_stamp
        self.noise_map = noise_map
        self.target_mask = target_mask
    
    def find_PSF(self, radius = 50):
        from decomprofile.tools_data.measure_tools import search_local_max, measure_FWHM
        init_PSF_locs_ = search_local_max(self.fov_image)
        init_PSF_locs, FWHMs, fluxs = [], [], []
        for i in range(len(init_PSF_locs_)):
            cut_image = cut_center_auto(self.fov_image, center = init_PSF_locs_[i], radius=radius)
            _fwhms = measure_FWHM(cut_image , radius = int(radius/5))
            if np.std(_fwhms)/np.mean(_fwhms) < 0.2 :
                init_PSF_locs.append(init_PSF_locs_[i])
                FWHMs.append(np.mean(_fwhms))
                fluxs.append(np.sum(cut_image))
        init_PSF_locs = np.array(init_PSF_locs)
        FWHMs = np.array(FWHMs)
        fluxs = np.array(fluxs)
        select_bool = (FWHMs<4.2)*(fluxs<5000)*(fluxs>200) #!!! A threshold to rule out the PSFs that are too board/bright/faint.
        PSF_locs = init_PSF_locs[select_bool]     
        for i in range(len(PSF_locs)):
            cut_image = cut_center_auto(self.fov_image, center = PSF_locs[i], kernel = 'center_gaussian', radius=radius)
            print('PSF location:', PSF_locs[i])
            print('id:', i, 'FWHMs:', np.round(measure_FWHM(cut_image , radius = int(radius/5)),3), 'flux:', round(np.sum(cut_image),1) )
            plt_fits(cut_image)
        select_idx = input('Input directly the a obj idx to mask, use space between each id:\n')
        select_idx = select_idx.split(" ")
        select_idx = [int(select_idx[i]) for i in range(len(select_idx)) if select_idx[i].isnumeric()]
        self.PSF_pos_list = [PSF_locs[i] for i in select_idx]
        self.PSF_lists = [cut_center_auto(self.fov_image, center = PSF_locs[i], kernel = 'center_gaussian', radius=radius) for i in select_idx]
        

# [] Test exptime in float or array        