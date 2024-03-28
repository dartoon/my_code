#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 15:17:43 2024

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle
from galight.tools.astro_tools import plt_fits
from galight.tools.astro_tools import plt_many_fits

filt =  'F150W'

idx = 9 #!!!

import sys
sys.path.insert(0, '../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info

info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

folder = '../../2022_JWST_QSOz6/NIRCam_data/*/bkg_removed/'   #!!!
from astropy.coordinates import SkyCoord
from astropy import units as u
# pos = SkyCoord('{0} {1}'.format(RA, Dec), unit=(u.hourangle, u.deg))
# target_pos = np.array([pos.ra.degree, pos.dec.degree])

jwst_all_filenames = glob.glob(folder+'*{0}*{1}*{2}*_rmbkg.fits'.format(target_id[:5],target_id[-4:], filt))  #For NIRCam
jwst_all_filenames.sort()

#%%
from galight.tools.cutout_tools import psf_clean
re_select = True
clean_up = True
save_folder = '../material/{0}/'.format(filt)


filename = jwst_all_filenames[0]
print("Select for", filename.split('/')[-1])
# Grab the JWST provided ERR map:

fitsFile = pyfits.open(filename)
fov_image = fitsFile[1].data # check the back grounp
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
    
flux_mjsr = header['PHOTMJSR']
pixscale = read_pixel_scale(header)
zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux

# fov_noise_map = fitsFile[2].data

wht = fitsFile[4].data # The WHT map
exp = fitsFile[0].header['EFFEXPTM']
gain_value = 2
exp_map = exp * wht/wht.max() / flux_mjsr * gain_value

from matplotlib.colors import LogNorm
data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
                          rm_bkglight = False, exptime = exp_map, if_plot=False, zp = zp)#, fov_noise_map = fov_noise_map)
data_process.find_PSF(radius = 60, user_option = True, if_filter=False, nearyby_obj_filter=False, FWHM_sort=True, 
                      norm = LogNorm(vmin = 0.001, vmax = 20))

remove_i = []
for i in range(len(data_process.PSF_list)):
    if np.sqrt(np.sum((data_process.PSF_pos_list[i] - data_process.target_pos)**2)) < 15:
        remove_i.append(i)
print("remove_i", remove_i)
data_process.PSF_list = [data_process.PSF_list[i] for i in range(len(data_process.PSF_list)) if i not in remove_i]
data_process.PSF_FWHM_list = [data_process.PSF_FWHM_list[i] for i in range(len(data_process.PSF_FWHM_list)) if i not in remove_i]
data_process.PSF_pos_list = [data_process.PSF_pos_list[i] for i in range(len(data_process.PSF_pos_list)) if i not in remove_i]
data_process.PSF_flux_list = [data_process.PSF_flux_list[i] for i in range(len(data_process.PSF_flux_list)) if i not in remove_i]
    
data_process.plot_overview()
from galight.tools.astro_tools import plt_many_fits


PSF_pos_list = data_process.PSF_pos_list
from astropy.wcs import WCS
wcs = WCS(header)
PSF_RA_DEC_list = []
for i in range(len(PSF_pos_list)):
    RA_DEC = wcs.all_pix2world(PSF_pos_list[i][0], PSF_pos_list[i][1],0) 
    PSF_RA_DEC_list.append( [float(RA_DEC[0]), float(RA_DEC[1])] )
data_process.PSF_RA_DEC_list = PSF_RA_DEC_list

PSF_list_clean = []
for i, psf in enumerate(data_process.PSF_list):
    print("Clean PSF", i)
    psf = data_process.PSF_list[i]
    psf = psf_clean(psf,if_plot=False, nsigma=3, npixels=45, ratio_to_replace=0.005,
                    if_print_fluxratio=True)
    PSF_list_clean.append(psf)
data_process.PSF_list_clean = PSF_list_clean
print("Before Clean")
plt_many_fits(data_process.PSF_list,norm = LogNorm(vmin = 0.001, vmax = 20))
print("After Clean")
plt_many_fits(data_process.PSF_list_clean,norm = LogNorm(vmin = 0.001, vmax = 20))
del data_process.fov_image
del data_process.exptime
pickle.dump(data_process, open(save_folder+'PSF_library_idx{0}_DPformat.pkl'.format(idx), 'wb'))
