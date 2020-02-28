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
sys.path.insert(0,'../../../py_tools/')
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result, transfer_obj_to_result
from mask_objects import detect_obj
from flux_profile import profiles_compare, flux_profile
from matplotlib.colors import LogNorm
import copy
import time
import pickle

#Setting the fitting condition:
deep_seed = True  #Set as True to put more seed and steps to fit.
pltshow = 1 #Note that setting plt.ion() in line27, the plot won't show anymore if running in terminal.
pix_scale = 0.04 
fixcenter = False
run_MCMC = False
zp= 28.

psf, QSO_img = pyfits.getdata('sim_ID_0/Drz_PSF.fits'),  pyfits.getdata('sim_ID_0/Drz_QSO_image.fits')

framesize = 60
ct = (len(QSO_img) - framesize)/2
QSO_img = QSO_img[ct:-ct,ct:-ct]

exptime = 625.0 * 8
stdd =  0.0088  #Measurement from empty retion.
QSO_std = (abs(QSO_img/exptime)+stdd**2)**0.5
#%%
#==============================================================================
# input the objects components and parameteres
#==============================================================================
objs, Q_index = detect_obj(QSO_img,pltshow = pltshow)
qso_info = objs[Q_index]
obj = [objs[i] for i in range(len(objs)) if i != Q_index]
fixed_source = []
kwargs_source_init = []
kwargs_source_sigma = []
kwargs_lower_source = []
kwargs_upper_source = []
fixed_source.append({})  
kwargs_source_init.append({'R_sersic': 0.3, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -0.5, 'center_y': -0.5})
kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': 0.5, 'center_y': 0.5})
#if len(obj) >= 1:
#    for i in range(len(obj)):
#        fixed_source.append({})  
#        kwargs_source_init.append({'R_sersic': obj[i][1] * pix_scale, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
#        kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
#        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': obj[i][1] * pix_scale/5, 'n_sersic': 0.3, 'center_x': -obj[i][0][0]*pix_scale-10, 'center_y': obj[i][0][1]*pix_scale-10})
#        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': -obj[i][0][0]*pix_scale+10, 'center_y': obj[i][0][1]*pix_scale+10})
source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]

#%%
#==============================================================================
# to fit and save the inference
#==============================================================================
tag = 'example'
source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_img, psf_ave=psf, psf_std = None,
                                                                  source_params=source_params, fixcenter=fixcenter,
                                                                  pix_sz = pix_scale, no_MCMC = (run_MCMC==False),
                                                                  QSO_std =QSO_std, tag=tag, deep_seed= deep_seed, pltshow=pltshow,
                                                                  corner_plot=False, flux_ratio_plot=True, dump_result=run_MCMC, pso_diag =True)

if pltshow == 0:
    plot_compare=False
    fits_plot =False
else:
    plot_compare=True
    fits_plot =True
    
result = transfer_to_result(data=QSO_img, pix_sz = pix_scale,
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, error_map=error_map,
                            zp=zp, fixcenter=fixcenter,ID='Example', tag=tag, plot_compare = plot_compare)


#%%
