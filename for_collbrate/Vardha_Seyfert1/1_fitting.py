#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 22:34:48 2019

@author: Dartoon

Fit the Seyfert as Point source + two Sersic (Bluge + Disk).
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'./fitting_tools/')
from fit_qso import fit_qso

ID = 'l106'

#Input the fitting ingredients:
agn_image = pyfits.getdata('{0}_cutout.fits'.format(ID))
agn_stdd = pyfits.getdata('{0}_stdd.fits'.format(ID))

##Activate this part if one need to cut to a smaller size such as 81:
#fit_frame_size = 81
#ct = (len(agn_image)-fit_frame_size)/2     # If want to cut to 81, agn_image[ct:-ct,ct:-ct]
#agn_image = agn_image[ct:-ct,ct:-ct]
#agn_stdd = agn_stdd[ct:-ct,ct:-ct]
print "The fitting image:"
plt.imshow(agn_image, norm = LogNorm(), origin='lower')
plt.show()

psf = pyfits.getdata('PSF0.fits')  # Input the PSF0 to fit.
##Activate this part if one want to use a smaller PSF size to make the fitting faster.
#psf_size = 81
#pct = (len(psf)-psf_size)/2     # If want to cut to 81, agn_image[ct:-ct,ct:-ct]
#psf = psf[pct:-pct,pct:-pct]
psf /= psf.sum()
print "The adopted PSF:"
plt.imshow(psf, norm = LogNorm(), origin='lower')
plt.show()

# =============================================================================
# Creat the QSO mask
# =============================================================================
from mask_objects import mask_obj
_ , _, deblend_sources = mask_obj(agn_image, snr=1, npixels=200, return_deblend = True)

print "deblend image to find the ID for the Objects for the mask:"
plt.imshow(deblend_sources, origin='lower',cmap=deblend_sources.cmap(random_state=12345))
plt.colorbar()
plt.show()

# Based on the deblend_sources, it is the ID 1 and 2 to block:
block_id = [1,2]
if block_id == []:
    QSO_msk = np.ones_like(agn_image)
else:
    for i in range(len(block_id)):
        if i ==0:
            mask = (np.array(deblend_sources)==block_id[i])
        else:
            mask += (np.array(deblend_sources)==block_id[i])
        
        QSO_msk = 1- mask
print "The QSO mask for the fitting:"
plt.imshow(QSO_msk, origin='lower')
plt.show()
# =============================================================================
# Setting up the parameter to fit
# =============================================================================
fixed_source = []
kwargs_source_init = []
kwargs_source_sigma = []
kwargs_lower_source = []
kwargs_upper_source = []      
fixed_source.append({})  
kwargs_source_init.append({'R_sersic': 2, 'n_sersic': 4., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.04, 'n_sersic': 0.5, 'center_x': -1, 'center_y': -1})
kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 16., 'n_sersic': 7,'center_x': 1, 'center_y': 1})
fixed_source.append({'n_sersic': 1.})  
kwargs_source_init.append({'R_sersic': 10, 'n_sersic': 1., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.04, 'center_x': -1, 'center_y': -1})
kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 16., 'center_x': 1, 'center_y': 1})
##Activate this part if one want to put another as Bar component:
#fixed_source.append({'n_sersic': 0.5})  
#kwargs_source_init.append({'R_sersic': 5, 'n_sersic': 0.5, 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
#kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
#kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.04, 'center_x': -1, 'center_y': -1})
#kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 16., 'center_x': 1, 'center_y': 1})
source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]

fixed_ps = []
kwargs_ps_init = []
kwargs_ps_sigma = []
kwargs_lower_ps = []
kwargs_upper_ps = []
fixed_ps.append({})
kwargs_ps_init.append({'ra_image': [0.], 'dec_image': [0.]})
kwargs_ps_sigma.append({'ra_image': [0.1], 'dec_image': [0.1]})
kwargs_lower_ps.append({'ra_image': [-1.], 'dec_image': [-1.]})
kwargs_upper_ps.append({'ra_image': [1.], 'dec_image': [1.]})
ps_param = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]         

name_save = 'fitting'
deep_seed = True
pix_sz = 0.04  # The pixel scale of the image.

source_result, ps_result, image_ps, image_host, error_map=fit_qso(agn_image, psf_ave=psf, ps_param = ps_param, 
                                                                  pix_sz = pix_sz, QSO_std = agn_stdd,  QSO_msk = QSO_msk,
                                                                  source_params=source_params, fixcenter=False, no_MCMC =False,
                                                                  tag=name_save, deep_seed= deep_seed, pltshow = 1, flux_ratio_plot=True)   
zp = 25.0985
print "Total Magnitude Result:"
print "AGN magnitude:", -2.5*np.log(image_ps.sum()) + zp
print "Buldge magnitude:", -2.5*np.log(image_host[0].sum()) + zp
print "Disk magnitude:", -2.5*np.log(image_host[1].sum()) + zp
#print "Bar magnitude:", -2.5*np.log(image_host[2].sum()) + zp
