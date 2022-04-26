#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 17:09:41 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
# import matplotlib.pyplot as plt

# import pandas as pd

# cata = pd.read_csv('TNG50_catalog/TNG50-1_091_Catalogue.csv')

# #Mass in log(M/M_sun), radius in kpc
# stellar_mass = cata['SubhaloMassType_stars']
# gas = cata['SubhaloMassType_gas']
# Reff = cata['SubhaloHalfmassRadType_stars']
# ID = cata['SubhaloID']
# ID_select = ID[(Reff>4)*(stellar_mass>10.5)] 

from galight.tools.plot_tools import plt_fits

# fits = pyfits.getdata('TNG50_img/shalo_091-544430_v3_photo.fits')
# fits[fits==99.] = 30
file = 'TNG50_img/shalo_091-544430_v3_photo.fits'
im = pyfits.open(file)
fits = im[0].data
header = im[0].header

zp=27
flux_hd = 10**(zp-fits)/0.4
from scipy import signal
from scipy.ndimage import zoom

z_s = 6.2
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
scale_relation = cosmo.angular_diameter_distance(z_s).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
scale = 0.1/scale_relation/ 0.03

# host_Reff = host_Reff_kpc/scale_relation   #In arcsec
flux = zoom(flux_hd, scale)

import pickle
_psfs, _FWHMs = pickle.load(open('../prep_use_HST_highRes/f356w_psfs.pkl','rb'))
psfs, FWHMs, fluxs = [], [], []
for i in range(len(_psfs)):
    psfs = psfs + _psfs[i]
psf = psfs[0]
flux = signal.fftconvolve(flux, psf, mode='full')
plt_fits( flux, norm = None  )