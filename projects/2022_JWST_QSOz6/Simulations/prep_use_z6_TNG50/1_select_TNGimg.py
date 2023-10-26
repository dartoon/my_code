#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 22:52:12 2022

@author: Dartoonw
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import glob
from scipy.ndimage import zoom
from galight.tools.astro_tools import plt_fits, plt_many_fits
from galight.tools.measure_tools import SB_profile, flux_profile

TNG_files = glob.glob("TNG50_z6/Idealized/013/*v0*fits")
TNG_files.sort()
TNG_ID = []
# TNG_files  = TNG_files[2:4]

zp = 27
pixscale = 0.03

band = ['JWST_NIRCAM.F150W', 'JWST_NIRCAM.F356W'][0]

# save_list = []
all_imgs = []
all_Reff_list = []
choose_imgs = []
choose_ids = []
# choose_ids = ['101491', '101499', '119454', '140982', '15', '219845', '272230', '321717', '561512', '579945']
choose_Reffs = []
for ID, name in enumerate(TNG_files):
    TNG_file = TNG_files[ID]
    im = pyfits.open(TNG_file)
    # fits = im['SUBARU_HSC.G'].data
    fits = im[band].data
    # _header = im['SUBARU_HSC.G'].header
    flux_hd = 10**(-0.4*(fits-zp))
    # plt_fits(flux_hd[450:-450, 450:-450])
    # flux_hd = flux_hd[450:-450, 450:-450]
    z_s = 6.2
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    scale_relation = cosmo.angular_diameter_distance(z_s).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/"
    scale = 0.1/scale_relation/ pixscale #Pixel scale is 100 pc, 0.1 kpc.
    # host_Reff = host_Reff_kpc/scale_relation   #In arcsec
    flux_zoom = zoom(flux_hd, scale)
    flux_zoom = flux_zoom/np.sum(flux_zoom) * 2000
    # plt_fits(flux_zoom)
    
    print(name)
    #Measrue the Reff:
    fluxs, rad, _ = flux_profile(flux_zoom, center=[len(flux_zoom)/2]*2 , radius=len(flux_zoom)/2, if_plot=True, #if_annuli=(True), 
                      fits_plot=(True), grids=50, x_gridspace=None)
    Reff_rad = rad[fluxs<fluxs[-1]/2][-1]
    Reff = Reff_rad / scale * 0.1  #Transfer to original fits's reff (0.1 kpc/pixel) and then to Kpc.

    SBs, _ = SB_profile(flux_zoom, center=[len(flux_zoom)/2]*2 , radius=len(flux_zoom)/2, if_plot=False, if_annuli=(True), 
                      fits_plot=(False), grids=50, x_gridspace=None)   
    all_imgs.append(flux_zoom)
    all_Reff_list.append(Reff)
    if TNG_file.split('-')[-1].split('_')[0] in choose_ids:
    # if Reff < 1.5 and Reff >0.5 and np.average(SBs[-4:])< 0.001:
        # print(name)
        # print(Reff_rad, Reff, '\n')
        # plt_fits(flux_zoom)
        # _ = flux_profile(flux_zoom, center=[len(flux_zoom)/2]*2 , radius=len(flux_zoom)/2, if_plot=True, #if_annuli=(True), 
        #                   fits_plot=(False), grids=50, x_gridspace=None)         
        # print(np.average(SBs[-4:]) )
        # save_list.append(name.split('/')[-1])
        choose_imgs.append(flux_zoom)
        choose_ids.append(TNG_file.split('-')[-1].split('_')[0])
        choose_Reffs.append(Reff)
        
        
TNG_ID = [TNG_files[i].split('-')[-1].split('_')[0] for i in range(len(TNG_files))]
plt_many_fits(all_imgs,labels = TNG_ID, prop = 'Reff(kpc)', texts = all_Reff_list)     
# plt_many_fits(choose_imgs,labels = choose_ids, prop = 'Reff(kpc)', texts = choose_Reffs, savename = 'select_'+band+'.pdf')     

