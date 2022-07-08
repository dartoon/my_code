#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 17:09:41 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

# import pandas as pd

# 	SubhaloFlag	SubhaloSFR	SubhaloSFRinHalfRad	SubhaloSFRinRad	SubhaloID	Snapshot	SubhaloHalfmassRadType_stars	SubhaloHalfmassRadType_gas	SubhaloHalfmassRadType_dm	SubhaloMassType_stars	SubhaloMassType_gas	SubhaloMassType_dm
# 544140	TRUE	5.144813	1.1004618	3.3023458	544140	91	9.414724	36.98669	68.459694	10.315633	10.892509	11.766261
# 544430	TRUE	0.42086697	0.0	0.047758806	544430	91	4.906846	78.38346	74.43424	10.5089855	10.413807	11.876655

# cata = pd.read_csv('TNG50_catalog/TNG50-1_091_Catalogue.csv')

# #Mass in log(M/M_sun), radius in kpc
# stellar_mass = cata['SubhaloMassType_stars']
# gas = cata['SubhaloMassType_gas']
# Reff = cata['SubhaloHalfmassRadType_stars']
# ID = cata['SubhaloID']
# ID_select = ID[(Reff>4)*(stellar_mass>10.5)] 

from galight.tools.plot_tools import plt_fits, plt_many_fits
from scipy import signal
from scipy.ndimage import zoom
from astropy.cosmology import FlatLambdaCDM
import pickle


# fits = pyfits.getdata('TNG50_img/shalo_091-544430_v3_photo.fits')
# fits[fits==99.] = 30
# file = 'TNG50_img/shalo_091-544140_v1_photo.fits'

import glob

files = glob.glob("TNG50_z6/Idealized/013/*7251*v?*fits")
files = glob.glob("TNG50_z6/Idealized/013/*v0*fits")
files.sort()
Reff_kpcs = []
for idx, file in enumerate(files): 
# for idx, file in enumerate(['TNG50_img/shalo_091-17_v0_photo.fits']):
    im = pyfits.open(file)
    # im.info()
    fits = im['JWST_NIRCAM.F356W'].data  #JWST_NIRCAM.F150W
    header = im['JWST_NIRCAM.F356W'].header
    
    # im['SUBARU_HSC.G'].header['WLPIVOT']*(1+6.2) #Compare to 3.56
    zp=22.5
    flux_hd = 10**((zp-fits)*0.4)
    z_s = 6.2
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    scale_relation = cosmo.angular_diameter_distance(z_s).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/"
    scale = 0.1/scale_relation/ 0.03 #Project 0.1 kpc to JWST pixel
    
    # host_Reff = host_Reff_kpc/scale_relation   #In arcsec
    flux_zoom = zoom(flux_hd, scale)
    
    _psfs, _FWHMs = pickle.load(open('../prep_use_HST_highRes/f356w_psfs.pkl','rb'))
    psfs, FWHMs, fluxs = [], [], []
    for i in range(len(_psfs)):
        psfs = psfs + _psfs[i]
    psf = psfs[1]
    psf = psf/psf.sum()
    flux = flux_zoom
    flux = signal.fftconvolve(flux_zoom, psf, mode='same')
    
    flux = flux/np.sum(flux) * 100
    flux_hd = flux_hd/np.sum(flux_hd)* 100
    # plt_fits(flux_hd)
    # plt_fits(flux)
    from matplotlib.colors import LogNorm
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 10))
    im1 = ax1.imshow(flux_hd, origin='lower', norm=LogNorm(vmax=flux_hd.max(), vmin = 1.e-5))
    ax1.set_title('HD image', fontsize=25)
    fig.colorbar(im1, ax=ax1, pad=0.01,  orientation="horizontal")
    # ax1.get_xaxis().set_visible(False)
    # ax1.get_yaxis().set_visible(False) 
    im2 = ax2.imshow(flux, origin='lower', norm=LogNorm(vmax=flux_hd.max(),  vmin = 1.e-5))
    ax2.set_title('Projected and convolved', fontsize=25)
    fig.colorbar(im2, ax=ax2, pad=0.01,  orientation="horizontal")
    # ax2.get_xaxis().set_visible(False)
    # ax2.get_yaxis().set_visible(False) 
    from galight.tools.measure_tools import SB_profile,flux_profile
    fluxs, rads, _ = flux_profile(flux_hd, center=[len(flux_hd)/2]*2 , radius=150, x_gridspace='log',
                                  if_plot=False, fits_plot=(False), grids = 40)
    try:
        Reff = rads[fluxs/100>0.5][0]
    except:
        Reff = rads[fluxs/fluxs[-1]>0.5][0]
    Reff_kpcs.append(Reff*0.1)
    print(file, idx, round(Reff*0.1,3), 'kpc:')
    plt.show()  
    
print("JWST PSF FWHM:", 100/scale * 4 / 1000, 'kpc')  #4 is the PSF typtical FWHM in pixels
    
    
    
    
    
    