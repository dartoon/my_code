#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 16:03:33 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.tools.plot_tools import plt_fits
Reines_t1 = np.loadtxt('2021_previous/Reines_2015_table_1.txt', dtype=str)

f = open("2021_previous/Reines_ID_RaDec.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n

mag_mis = []
IDs = glob.glob('s21a/*')

images = []
IDs.sort()
for ID in IDs:
    ID = ID.split('/')[1]
    idx = np.where(Reines_t1[:,2] == ID)[0][0]
    Reines_iMag = float(Reines_t1[idx, 5])
    z = float(Reines_t1[idx, 4])
    
    ra = ID[1:10]
    ra = ra[:2] + ':' + ra[2:4] + ':' +  ra[4:]
    dec = ID[10:]
    dec = dec[:3] + ':' + dec[3:5] + ':' +  dec[5:]
    
    fitsname = glob.glob('./s21a/{0}/*cutout*I*.fits'.format(ID))
    fitsFile = pyfits.open(fitsname[0])
    fov_image= fitsFile[1].data
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    err_data= fitsFile[3].data ** 0.5
    
    file_header0 = fitsFile[0].header
    FLUXMAG0 = file_header0['FLUXMAG0']
    zp =  2.5 * np.log10(FLUXMAG0)
    
    psfname = glob.glob('./s21a/{0}/*psf*I*.fits'.format(ID))
    PSF = pyfits.getdata(psfname[0])
    

    data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [ra, dec],
                               pos_type = 'wcs', header = header,
                               rm_bkglight = True, if_plot=False, zp = zp)
    data_process.generate_target_materials(radius=180, create_mask = False, nsigma=2.8,
                                           radius_list = [80, 100, 120, 140, 160, 180], 
                                          exp_sz= 1.2, npixels = 15, if_plot=False)
    data_process.PSF_list = [PSF]
    data_process.checkout() #Check if all the materials is known.

    fit_sepc = FittingSpecify(data_process)
    fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor=3)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
    # plt_fits(data_process.target_stamp)
    images.append(data_process.target_stamp)
    # fit_sepc.plot_fitting_sets()
    # fit_sepc.build_fitting_seq()

    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    dl = cosmo.luminosity_distance(z).value  # Mpc
    obs_mag =  -2.5*np.log10(np.sum(data_process.target_stamp)) + 27
    abs_mag = obs_mag - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
    print("mag miss:", abs_mag - Reines_iMag)
    mag_mis.append(abs_mag - Reines_iMag)
    
#%%
from galight.tools.astro_tools import plt_many_fits
plt_many_fits(images, labels = [IDs[i][5:] for i in range(len(IDs))], label_size=13)

    
    