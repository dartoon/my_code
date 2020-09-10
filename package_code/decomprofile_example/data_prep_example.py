#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 16:01:17 2020

@author: Xuheng Ding

You can skip this step if the QSO stamp, noise level and the PSF is ready.
"""
#photutils in version 0.7.2
#astropy in version astropy-4.0.1


import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from matplotlib.colors import LogNorm

from decomprofile.tools_data.astro_tools import plt_fits
#Load data and plot:
fitsFile = pyfits.open('../example_data/HST/QSO/1104_final_drz.fits')
img = fitsFile[1].data # check the back grounp

plt_fits(img, figsize=(15,15))

#%%
import decomprofile.tools_data.astro_tools as astro_tools
print(astro_tools.read_fits_exp(fitsFile), astro_tools.read_pixel_scale(fitsFile,frame=1))  #Read the exposure time and pixel scale.

from decomprofile.tools_data.measure_tools import measure_bkg
bkglight = measure_bkg(img, if_plot=True)
img = img - bkglight   #remove the bkglight

#%%
from decomprofile.tools_data.cutout_tools import cut_center_auto
QSO_loc = [1135, 648]  # The postion of the QSO in the frame
QSO_img, QSO_center_pos = cut_center_auto(image=img, center= QSO_loc, 
                                          kernel = 'center_gaussian', radius=60,
                                          return_center=True, if_plot=True)
plt_fits(QSO_img, colorbar = True)

#%%Creat the SB profile of the QSO:
from decomprofile.tools_data.measure_tools import SB_profile, esti_bgkstd
# r_SB, r_grids = SB_profile(QSO_img, center = [(len(QSO_img)-1)/2]*2 , radius=20,
#                            grids=50, x_gridspace='log',if_annuli=True, if_plot=True, fits_plot = True)
std = esti_bgkstd(QSO_img, if_plot=True)

#%%Test the way to creat the mask for the QSO:
from decomprofile.tools_data.measure_tools import detect_obj, mask_obj
apertures = detect_obj(QSO_img, if_plot=True)

select_mask_idx = [0,1,3]

apertures = [apertures[i] for i in select_mask_idx]
mask_list = mask_obj(QSO_img, apertures, if_plot=False)
mask = np.ones_like(QSO_img)
for i in range(len(mask_list)):
    mask *= mask_list[i]
plt_fits(mask, colorbar = False)

#%%Auto find the PSF in the frames
from decomprofile.tools_data.measure_tools import search_local_max, measure_FWHM
from decomprofile.tools_data.cutout_tools import cutout
init_PSF_locs = search_local_max(img)
init_PSF_locs = np.array(init_PSF_locs)
FWHMs, fluxs = [], []
for i in range(len(init_PSF_locs)):
    cut_image = cut_center_auto(img, center = init_PSF_locs[i], radius=60)
    FWHMs.append(np.mean(measure_FWHM(cut_image , radius = 10)))
    fluxs.append(np.sum(cut_image))
FWHMs = np.array(FWHMs)
fluxs = np.array(fluxs)
select_bool = (FWHMs<4.2)*(fluxs<5000)*(fluxs>200) # A threshold to rule out the PSFs that are too board/bright/faint.
PSF_locs = init_PSF_locs[select_bool]

#%%
for i in range(len(PSF_locs)):
    cut_image = cut_center_auto(img, center = PSF_locs[i], kernel = 'center_gaussian', radius=60)
    print('PSF location:', PSF_locs[i])
    print('id:', i, 'FWHMs:', np.round(measure_FWHM(cut_image , radius = 10),3), 'flux:', round(np.sum(cut_image),) )
    plt_fits(cut_image)
    
select = [3, 4, 5, 10]
PSF_locs_final = [PSF_locs[i] for i in select]
PSF_lists = [cut_center_auto(img, center = PSF_locs[i], kernel = 'center_gaussian', radius=50) for i in select]


from decomprofile.tools_data.measure_tools import profiles_compare
profiles_compare([QSO_img] + PSF_lists, x_gridspace= 'log', norm_pix = 5, if_annuli=True, y_log = False,
                 prf_name_list = (['QSO'] + ['PSF{0}'.format(i) for i in range(len(PSF_lists))]))

from decomprofile.tools_data.cutout_tools import plot_overview
plot_overview(img, center_QSO= QSO_center_pos, c_psf_list=PSF_locs_final)
pyfits.PrimaryHDU(QSO_img).writeto('QSO_image.fits', overwrite=True)
for i in range(len(PSF_locs_final)):
    PSF_cut = cutout(img, center = PSF_locs_final[i], radius=50)
    pyfits.PrimaryHDU(PSF_cut).writeto('PSF{0}.fits'.format(i),overwrite=True)
