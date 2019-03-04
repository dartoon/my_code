#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 22:34:48 2019

@author: Dartoon

Fit the double as two Point source + a Sersic.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools/')
from fit_qso import fit_qso

ID = 'l106'

agn_image = pyfits.getdata('{0}_cutout.fits'.format(ID))
agn_stdd = pyfits.getdata('{0}_stdd.fits'.format(ID))
psf = pyfits.getdata('PSF0.fits')
psf /= psf.sum()

##If need to cut 
#fit_frame_size = 301
#ct = (len(agn_image)-fit_frame_size)/2     # If want to cut to 61, QSO_im[ct:-ct,ct:-ct]
#agn_image = agn_image[ct:-ct,ct:-ct]
#agn_rms = agn_rms[ct:-ct,ct:-ct]

plt.imshow(agn_image, norm = LogNorm(), origin='lower')
plt.show()

pix_sz = 0.04

#from rebin import block, expend_grid
#agn_image = block(expend_grid(agn_image,3), shape=[201,201],factor=4)
#plt.imshow(agn_image, origin='lower', norm=LogNorm())
#plt.colorbar()
#plt.show()
#from photutils import make_source_mask
#mask = make_source_mask(agn_image, snr=1, npixels=5, dilate_size=11)
#plt.imshow(agn_image * (1-mask*1), norm=LogNorm(), origin='low')
#plt.show()
#stdd = np.std(agn_image[mask ==False])  # 
#print stdd
#agn_expmap= pyfits.getdata('{0}_exp_map.fits'.format(ID))
#agn_expmap = block(expend_grid(agn_expmap,3), shape=[201,201],factor=4)
#agn_stdd = (abs(agn_image/agn_expmap)+stdd**2)**0.5
#psf = block(psf[3:,3:], shape=[37,37],factor=4)
#pix_sz = pix_sz*4

# =============================================================================
# Creat the mask
# =============================================================================
from mask_objects import mask_obj
_ , _, deblend_sources = mask_obj(agn_image, snr=1, npixels=200, return_deblend = True)

plt.imshow(deblend_sources, origin='lower',cmap=deblend_sources.cmap(random_state=12345))
plt.colorbar()
plt.show()

QSO_msk = np.ones_like(agn_image)
QSO_msk = 1- ((np.array(deblend_sources)==2) + (np.array(deblend_sources)==1))
#QSO_msk = 1- (np.array(deblend_sources)==2)
print "The QSO mask"
plt.imshow(QSO_msk, origin='lower')
plt.show()

# =============================================================================
# Setting up the parameter to fit
# =============================================================================
#fixed_source = []
#kwargs_source_init = []
#kwargs_source_sigma = []
#kwargs_lower_source = []
#kwargs_upper_source = []      
#fixed_source.append({})  
#kwargs_source_init.append({'R_sersic': 4, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
#kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
#kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.04, 'n_sersic': 0.3, 'center_x': -1, 'center_y': -1})
#kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 16., 'n_sersic': 7., 'center_x': 1, 'center_y': 1})

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
fixed_source.append({'n_sersic': 0.5})  
kwargs_source_init.append({'R_sersic': 5, 'n_sersic': 0.5, 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.04, 'center_x': -1, 'center_y': -1})
kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 16., 'center_x': 1, 'center_y': 1})


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

source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
ps_param = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]         
   
source_result, ps_result, image_ps, image_host, error_map=fit_qso(agn_image, psf_ave=psf, ps_param = ps_param, 
                                                                  pix_sz = pix_sz, QSO_std = agn_stdd,  QSO_msk = QSO_msk,
                                                                  source_params=source_params, fixcenter=False, no_MCMC =True,
                                                                  tag='fitting_include_bar', deep_seed= True, pltshow = 1, flux_ratio_plot=True)   
zp = 25.0985
print "Total Magnitude Result:"
print "AGN magnitude:", -2.5*np.log(image_ps.sum()) + zp
print "Buldge magnitude:", -2.5*np.log(image_host[0].sum()) + zp
print "Disk magnitude:", -2.5*np.log(image_host[1].sum()) + zp
print "Bar magnitude:", -2.5*np.log(image_host[2].sum()) + zp
