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
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result, transfer_obj_to_result
from mask_objects import detect_obj
from flux_profile import profiles_compare, flux_profile
from matplotlib.colors import LogNorm
import copy
import time
import pickle

#Setting the fitting condition:
deep_seed = False  #Set as True to put more seed and steps to fit.
pltshow = 1 #Note that setting plt.ion() in line27, the plot won't show anymore if running in terminal.
pix_scale = 0.0642 
fixcenter = False
run_MCMC = True
zp=26.4524

psf, QSO_img, QSO_std = pyfits.getdata('PSF.fits'),  pyfits.getdata('QSO_im.fits'),  pyfits.getdata('QSO_err.fits')
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
# =============================================================================
# Creat the QSO mask
# =============================================================================
from mask_objects import mask_obj
_ , _, deblend_sources = mask_obj(QSO_img, snr=1.3, npixels=20, return_deblend = True)

print("deblend image to find the ID for the Objects for the mask:")
plt.imshow(deblend_sources, origin='lower',cmap=deblend_sources.cmap(random_state=12345))
plt.colorbar()
if pltshow == 0:
    plt.close()
else:
    plt.show()
block_id = [1]
if block_id == []:
    QSO_msk = np.ones_like(QSO_img)
else:
    for i in range(len(block_id)):
        if i ==0:
            mask = (np.array(deblend_sources)==block_id[i])
        else:
            mask += (np.array(deblend_sources)==block_id[i])
        
        QSO_msk = 1- mask
print("The QSO mask for the fitting:")
plt.imshow(QSO_msk, origin='lower')
if pltshow == 0:
    plt.close()
else:
    plt.show()
#%%
#==============================================================================
# to fit and save the inference
#==============================================================================
tag = 'example'
source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_img, psf_ave=psf, psf_std = None,
                                                                  source_params=source_params, QSO_msk = QSO_msk, fixcenter=fixcenter,
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
                            zp=zp, fixcenter=fixcenter,ID='Example', QSO_msk = QSO_msk, tag=tag, plot_compare = plot_compare)


#%%
