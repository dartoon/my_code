#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 16:35:06 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

#Load a fit_run
import glob, pickle
filename = '/Users/Dartoon/Astro/Projects/my_code/projects/2022_JWST_QSOz6/model_z6_data_id6/stage3_all/fit_material/fit_run_idx6_F356W_CombPsfsNO_2_0.pkl'
fit_run = pickle.load(open(filename,'rb'))
fit_spec = fit_run.fitting_specify_class


#%%
import os
import shutil

path = "galfit_convert/"
# Check whether the specified path exists or not
isExist = os.path.exists(path)
if isExist:
    shutil.rmtree(path)
os.makedirs(path)
#%%
data = fit_spec.kwargs_data['image_data']
pyfits.PrimaryHDU(data).writeto(path+'data.fits',overwrite=True)

noise = fit_spec.kwargs_data['noise_map']
pyfits.PrimaryHDU(noise).writeto(path+'noise.fits',overwrite=True)

psf = fit_spec.kwargs_psf['kernel_point_source']
pyfits.PrimaryHDU(psf).writeto(path+'psf.fits',overwrite=True)


#%%
f = open("string0.txt","r")
string_0 = f.read()
string_0 = string_0.replace('var_datafits', 'data.fits')
string_0 = string_0.replace('var_noisefits','noise.fits')
string_0 = string_0.replace('var_psffits', 'psf.fits')
string_0 = string_0.replace('var_size', str(fit_spec.numPix))
string_0 = string_0.replace('var_conv_size', str(len(fit_spec.psf_class.kernel_pixel)))
string_0 = string_0.replace('var_zp', str(fit_spec.zp))
string_0 = string_0.replace('var_resolution', str(fit_spec.deltaPix))
string_0 = string_0 + '\n'

#%%
from photutils.aperture import aperture_photometry
apertures = fit_spec.data_process_class.apertures
# i = 0
string_gal = ''
for i in range(len(apertures)):
    aper = apertures[i]
    f = open("string1_gal.txt","r")
    string_gal_i = f.read()
    var_objno = i + 1
    var_posx = aper.positions[0]
    var_posy = aper.positions[1]
    flux_in_ap = aperture_photometry(fit_spec.kwargs_data['image_data'], aper)['aperture_sum'].value[0]
    var_ini_mag = -2.5*np.log10(flux_in_ap) + fit_spec.zp
    var_ini_reff = aper.a
    var_ini_n = 2
    var_ini_q = aper.b / aper.a
    var_ini_theta = aper.theta /np.pi  * 180
    string_gal_i = '\n' + string_gal_i
    string_gal_i = string_gal_i.replace('var_objno', str(var_objno))
    string_gal_i = string_gal_i.replace('var_posx', str(var_posx))
    string_gal_i = string_gal_i.replace('var_posy', str(var_posy))
    string_gal_i = string_gal_i.replace('var_ini_mag', str(var_ini_mag))
    string_gal_i = string_gal_i.replace('var_ini_reff', str(var_ini_reff))
    string_gal_i = string_gal_i.replace('var_ini_n', str(var_ini_n))
    string_gal_i = string_gal_i.replace('var_ini_q', str(var_ini_q))
    string_gal_i = string_gal_i.replace('var_ini_theta', str(var_ini_theta))
    string_gal_i = string_gal_i + '\n'
    string_gal = string_gal+string_gal_i
#%%

string_ps = ''
if 'point_source_model' in fit_spec.kwargs_params.keys():
    for i in range(len(fit_spec.kwargs_params['point_source_model'][0])):
        f = open("string1_ps.txt","r")
        string_ps_i = f.read()
        var_objno = i + len(apertures) + 1
        var_posx = -fit_spec.kwargs_params['point_source_model'][0][i]['ra_image'][0] / fit_spec.deltaPix + fit_spec.numPix/2
        var_posy = fit_spec.kwargs_params['point_source_model'][0][i]['dec_image'][0] / fit_spec.deltaPix + fit_spec.numPix/2
        var_ini_mag = -2.5*np.log10(np.max(fit_spec.kwargs_data['image_data'])) + fit_spec.zp
        
        string_ps_i = '\n' + string_ps_i
        string_ps_i = string_ps_i.replace('var_objno', str(var_objno))
        string_ps_i = string_ps_i.replace('var_posx', str(var_posx))
        string_ps_i = string_ps_i.replace('var_posy', str(var_posy))
        string_ps_i = string_ps_i.replace('var_ini_mag', str(var_ini_mag))
        string_ps_i = string_ps_i + '\n'
        string_ps = string_ps +  string_ps_i
#%%

string_final = string_0 + string_gal + string_ps
write_file = open(path+'galfit.feedme','w') 
write_file.write(string_final)
write_file.close()
