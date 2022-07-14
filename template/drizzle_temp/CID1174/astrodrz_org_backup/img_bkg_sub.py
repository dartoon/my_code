#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 14:04:02 2018

@author: Dartoon

Subtract the bkg for each individual image
"""
import numpy as np
import sys
sys.path.insert(0,'../../../../py_tools')
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import glob

from astropy.visualization import SqrtStretch
from astropy.stats import SigmaClip
from photutils import Background2D, SExtractorBackground  
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import make_source_mask

files = glob.glob('*flt.fits')
for i in range(len(files)):
    fitsFile = pyfits.open(files[i], mode='update')
    img = fitsFile[1].data
    norm = ImageNormalize(stretch=SqrtStretch())         
    sigma_clip = SigmaClip(sigma=3., maxiters=10)
    bkg_estimator = SExtractorBackground()
    mask_0 = make_source_mask(img, nsigma=2, npixels=5, dilate_size=11)
    #mask_1 = (np.isnan(img))
    mask = mask_0 #+ mask_1
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    ax.imshow(mask,origin='lower')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False) 
    plt.show()
    bkg = Background2D(img, (50, 50), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                       mask=mask)
#    plt.imshow(img, norm=LogNorm(), origin='lower') 
#    #bkg.plot_meshes(outlines=True, color='#1f77b4')
#    plt.show()       
    back = bkg.background#* ~mask_1
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    ax.imshow(back, origin='lower', cmap='Greys_r')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)     
    plt.close()
    
    img -= back
    fitsFile[1].data = img
    fitsFile.flush()
    
print("total drz:", len(files))
    
             