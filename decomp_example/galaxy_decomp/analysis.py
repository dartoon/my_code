#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:58:39 2019

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
#from gen_fit_id import gen_fit_id
from photutils import make_source_mask
import os
#plt.ion()

import sys
sys.path.insert(0,'../../py_tools/')
from fit_qso import fit_galaxy
from transfer_to_result import transfer_to_result, transfer_obj_to_result
from mask_objects import detect_obj
from flux_profile import profiles_compare, flux_profile, galaxy_total_compare
from matplotlib.colors import LogNorm
import copy
import time
import pickle

#Setting the fitting condition:
deep_seed = False  #Set as True to put more seed and steps to fit.
pltshow = 1 #Note that setting plt.ion() in line27, the plot won't show anymore if running in terminal.
pix_scale = 0.168
run_MCMC = True
zp = 27.0

psf, galaxy_img, galaxy_std = pyfits.getdata('PSF.fits'),  pyfits.getdata('galaxy_im.fits'),  pyfits.getdata('galaxy_err.fits')
#%%
#==============================================================================
# input the objects components and parameteres
#==============================================================================
objs, Q_index = detect_obj(galaxy_img,pltshow = pltshow, snr=1.2)
galaxy_info = objs[Q_index]
obj = [objs[i] for i in range(len(objs)) if i != Q_index]
fixed_source = []
kwargs_source_init = []
kwargs_source_sigma = []
kwargs_lower_source = []
kwargs_upper_source = []
fixed_source.append({})  
kwargs_source_init.append({'R_sersic': 1, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -galaxy_info[0][0]*pix_scale, 'center_y': galaxy_info[0][1]*pix_scale})
kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -0.5, 'center_y': -0.5})
kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10., 'n_sersic': 7., 'center_x': 0.5, 'center_y': 0.5})
if len(obj) >= 1:
    for i in range(len(obj)):
        fixed_source.append({})  
        kwargs_source_init.append({'R_sersic': 1, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
        kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -obj[i][0][0]*pix_scale-5*pix_scale, 'center_y': obj[i][0][1]*pix_scale-5*pix_scale})
        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': -obj[i][0][0]*pix_scale+5*pix_scale, 'center_y': obj[i][0][1]*pix_scale+5*pix_scale})

source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
#%%
#==============================================================================
# to fit and save the inference
#==============================================================================
tag = 'example'
source_result, image_host, error_map, reduced_Chisq=fit_galaxy(galaxy_img, psf_ave=psf, psf_std = None,
                                                      source_params=source_params, galaxy_msk = None,
                                                      pix_sz = pix_scale, no_MCMC = (run_MCMC==False),
                                                      galaxy_std =galaxy_std, tag=tag, deep_seed= deep_seed,
                                                      pltshow=pltshow, return_Chisq=True,  corner_plot=False, 
                                                      dump_result=run_MCMC, flux_corner_plot = True, pso_diag=True)
if pltshow == 0:
    plot_compare=False
    fits_plot =False
else:
    plot_compare=True
    fits_plot =True

objs_img = np.zeros_like(image_host[0])
if len(source_result)>1:
    for i in range(1,len(image_host)):
        objs_img += image_host[i]    
flux_list = [galaxy_img, image_host[0]*0, image_host[0]+objs_img, error_map]    
result = transfer_obj_to_result(source_result[0],image_host[0], zp)
fig = galaxy_total_compare(label_list = ['data', 'model', 'normalized residual'], flux_list = flux_list, target_ID = 'example', pix_sz=pix_scale,
          zp = zp, plot_compare=plot_compare, msk_image = np.ones_like(image_host[0]))


#%%
