#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 17:28:30 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle
import glob
from astropy.cosmology import FlatLambdaCDM


seed = 5
# target_info = pickle.load(open('target_info.pkl','rb'))

# target_info = pickle.load(open('target_info.pkl','rb'))
res_files = glob.glob('sim_results/qsoID*filt_f35*seed{0}*result.pkl'.format(seed))
res = pickle.load(open(res_files[0],'rb'))

# file_fit = glob.glob('sim_results/qsoID0_filt_f356w_seed{0}diff_psf_fit.pkl'.format(seed))[0]
file_fit = glob.glob('sim_results/qsoID0_filt_f356w_seed{0}same_psf_fit.pkl'.format(seed))[0]
fit_run = pickle.load(open(file_fit,'rb'))
fit_run.plot_final_qso_fit()
print("Host ratio:", res['true_host_flux_ratio'])
print( "Inf flux VS True flux:", round(res['inferred_host_flux_same_psf'], 2), round( res['true_host_flux'], 2) )

Reff_pixel = fit_run.final_result_galaxy[0]['R_sersic'] / fit_run.fitting_specify_class.deltaPix

z_s = 6.2
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
scale_relation = cosmo.angular_diameter_distance(z_s).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/"
scale = 0.1/scale_relation/ 0.03 #Project 0.1 kpc to JWST pixel
   
inf_Reff_kpc = Reff_pixel *  100/scale/ 1000

if seed == 3:
    file = 'TNG50_img/shalo_091-540258_v0_photo.fits'
elif seed == 4:
    file = 'TNG50_img/shalo_091-542669_v0_photo.fits'
elif seed == 5:
    file = 'TNG50_img/shalo_091-572599_v0_photo.fits'

im = pyfits.open(file)
fits = im['SUBARU_HSC.G'].data
zp=27
flux_hd = 10**((zp-fits)*0.4)
from galight.tools.measure_tools import SB_profile,flux_profile
fluxs, rads, _ = flux_profile(flux_hd, center=[len(flux_hd)/2]*2 , radius=150, x_gridspace='log',
                             if_plot=False, fits_plot=(False), grids = 40)
True_Reff_kpc = rads[fluxs/np.sum(flux_hd)>0.5][0] * 0.1
print("True Reff (kpc) VS, inferred Reff", round(True_Reff_kpc,2), round(inf_Reff_kpc,2))