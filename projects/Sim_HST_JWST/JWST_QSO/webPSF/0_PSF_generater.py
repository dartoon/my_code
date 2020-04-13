#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 16:11:45 2020

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
from astropy.cosmology import FlatLambdaCDM
import glob
import webbpsf
#sys.path.insert(0, '../share_tools')

#%%Example of using grid:
#https://github.com/spacetelescope/webbpsf/blob/master/notebooks/Gridded_PSF_Library.ipynb
## Create a 3x3 grid of PSFs for NIRCam
#nrc = webbpsf.NIRCam()
#nrc.filter = "F444W"
#nrc_grid = nrc.psf_grid(num_psfs=9, all_detectors=False, oversample=4)
#webbpsf.gridded_library.display_psf_grid(nrc_grid)
#print(nrc_grid.meta)

#%%Set up data basic information
filt_l = ['F444W', 'F356W', 'F200W', 'F150W']
filt  = filt_l[3]
oversample = 4
position_list= [(1024.0, 1024.0),
        (0.0, 0.0),
      (0.0, 1024.0),
      (0.0, 2047.0),
      (1024.0, 0.0),
      (1024.0, 2047.0),
      (2047.0, 0.0),
      (2047.0, 1024.0),
      (2047.0, 2047.0)]
nc = webbpsf.NIRCam()
nc.filter =  filt
psf_list = []
for i in range(len(position_list)):
    print("Simulate PSF", i)
    nc.detector_position = position_list[i]
    psf_list.append(nc.calc_psf(oversample=oversample))     # returns an astropy.io.fits.HDUlist containing PSF and header
#plt.imshow((psf_list[0][0].data-psf_list[4][0].data), origin='lower')#,cmap='bwr')
#plt.colorbar()
#plt.show()    
#%%Create the saving file folder
#import os
#PSF_folder_name = 'highres_PSF_' + filt
#if os.path.exists(PSF_folder_name)==True:
#    ifdel = input("Simulation of "+ PSF_folder_name + " exist, delete the old folder?(Y/N):\n")
#    if ifdel == 'Y':
#        import shutil 
#        shutil.rmtree(PSF_folder_name)
#    elif ifdel == 'N':
#        print("Simulation stop!")
#        from sys import exit
#        exit
#    else:
#        raise ValueError("Input is not right, stop simulation.")        
#os.mkdir(PSF_folder_name)
#for i in range(len(psf_list)):
#    pyfits.PrimaryHDU(psf_list[i][0].data,header=psf_list[i][0].header).writeto(PSF_folder_name+'/'+'PSF_id{0}.fits'.format(i))
##    
