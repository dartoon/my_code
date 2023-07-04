#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 17:43:35 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob
import sys
import pickle
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.tools.measure_tools import measure_bkg
from astropy.wcs import WCS
from galight.tools.cutout_tools import common_data_class_aperture
from galight.tools.plot_tools import plot_data_apertures_point
from galight.tools.cutout_tools import cutout
import warnings
warnings.filterwarnings("ignore")
#%%
import sys
sys.path.insert(0,'../model_z6_data_id0/')

from target_info import target_info
from galight.tools.astro_tools import plt_fits, plt_many_fits

import copy, matplotlib
from matplotlib.colors import LogNorm

cmap = 'gist_heat'
my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap)) # copy the default cmap
my_cmap.set_bad('black')

data_type = 'all' 
filt = 'F356W'

file_NO = 0
idx = 1
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
folder = '../NIRCam_data/Nov14/bkg_removed/' 
jwst_all_filenames = glob.glob(folder+'/*{0}*.fits'.format(filt))
jwst_all_filenames.sort()
file = jwst_all_filenames[file_NO]

_fitsFile = pyfits.open(file)
fov_image = _fitsFile[1].data # check the back grounp
header = _fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
flux_mjsr = header['PHOTMJSR']
pixscale = read_pixel_scale(header)
zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
wht = _fitsFile[4].data # The WHT map
exp = _fitsFile[0].header['EFFEXPTM']
print("Exp time:", exp)
if _fitsFile[0].header['CHANNEL'] == 'LONG':
    gain_value = 2
    expsize = 1
    exppix = 1
else:
    gain_value = 1.8
    expsize = 1
    exppix = 2
exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], #The final cut center is dete
                               pos_type = 'wcs', header = header,rm_bkglight = False, 
                               if_plot=False, zp = zp, exptime= exp_map, 
                               fov_noise_map = None)
data_process.fov_noise_map  = data_process.fov_image
data_process.generate_target_materials(radius=800 * expsize, create_mask = False, nsigma=1.5, 
                                        cut_kernel = None, skip = True)
fig, ax = plt.subplots(figsize=(12,12))  #!!!
vmin = 1.e-3
if filt == 'F150W':
    vmax = data_process.target_stamp.max()/3 
else:
    vmax = data_process.target_stamp.max()
    
vmax = 9.4
plt.imshow(data_process.target_stamp, origin='lower', cmap=my_cmap, norm=LogNorm(vmin=vmin, vmax=vmax))#, vmin=vmin, vmax=vmax)
ax.set_xticks([])
ax.set_yticks([])
frame_size = len(data_process.target_stamp)
# plt.savefig('/Users/Dartoon/Downloads/J2236_fov.pdf')
plt.show()

#%%

data_process.generate_target_materials(radius=60 * expsize, create_mask = False, nsigma=1.5, 
                                        cut_kernel = None, skip = True)
fig, ax = plt.subplots(figsize=(4,4))  #!!!
vmax = data_process.target_stamp.max()
plt.imshow(data_process.target_stamp, origin='lower', cmap=my_cmap, norm=LogNorm(vmin=vmin, vmax=vmax))#, vmin=vmin, vmax=vmax)
ax.set_xticks([])
ax.set_yticks([])
frame_size = len(data_process.target_stamp)
# plt.savefig('/Users/Dartoon/Downloads/J2236_data.pdf')
plt.show()
#%%
run_folder = 'stage3_all/' #!!!
fit_files = glob.glob(run_folder+'*fit_material*/fit_run_fixn1__idx{0}_{1}_*.pkl'.format(idx, filt))#+\
# fit_files = glob.glob(run_folder+'*fit_material*/fit_run_withcentralMask_idx{0}_{1}_FOV*.pkl'.format(idx, filt))#+\
# fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
fit_files.sort()
fit_run_list=[]
for i in range(len(fit_files)):
    fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
sort_Chisq = chisqs.argsort()  
fit_run = fit_run_list[sort_Chisq[0]]
fig, ax = plt.subplots(figsize=(4,4))  #!!!
plt.imshow(fit_run.flux_2d_out['data-point source'], origin='lower', cmap=my_cmap, norm=LogNorm(vmin=vmin, vmax=vmax))#, vmin=vmin, vmax=vmax)
ax.set_xticks([])
ax.set_yticks([])
# plt.savefig('/Users/Dartoon/Downloads/J2236_data-qso.pdf')
#%%
plt.show()
fig, ax = plt.subplots(figsize=(4,4))  #!!!
plt.imshow(fit_run.image_ps_list[0], origin='lower', cmap=my_cmap, norm=LogNorm(vmin=0.0024, vmax=vmax/3))#, vmin=vmin, vmax=vmax)
ax.set_xticks([])
ax.set_yticks([])
plt.savefig('/Users/Dartoon/Downloads/J2236_qso.pdf')
plt.show()


