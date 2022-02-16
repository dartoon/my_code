#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 13:29:40 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

filt = 'f356w'
folder = 'JWST_CEERS/'
file = 'ceers5_{filt}_i2d.fits'.format(filt=filt)

im = pyfits.open(folder+file)
data_sb = im[1].data
header = im[1].header

data_sb = data_sb[:2000,:2000]
data = data_sb/header['PHOTMJSR']

value_unit = header['BUNIT']
print('For flux value in unit of MJy/sr, the zp is 31.4.') #https://en.wikipedia.org/wiki/Jansky#AB_magnitude

# print(header['BUNIT'])
print("Conversion factor from {units} to DN/S for filter {f}:".format(units=header['BUNIT'], f=filt), 
      header['PHOTMJSR'])

header0 = im[0].header
img_filter = header0['FILTER']
img_cam = header0['APERNAME'] #In JDAT'simulation it is 'DETECTOR'
exptime = header0['TEXPTIME']
from galight.tools.astro_tools import plt_fits
plt_fits(data)

# from galight.tools.measure_tools import detect_obj
# apr = detect_obj(data, if_plot=True, npixels=40)
#%% Model using galight
from matplotlib.colors import LogNorm
from galight.tools.plot_tools import my_cmap
#Start to use galight
from galight.data_process import DataProcess
zp = 31.4 - 2.5*np.log10(header['PHOTMJSR'])
data_process = DataProcess(fov_image = data, target_pos = [1170., 940.], pos_type = 'pixel', header = header,
                          rm_bkglight = False, exptime = np.ones_like(data)*exptime, if_plot=False, zp = zp)
data_process.generate_target_materials(radius=65, create_mask = False, nsigma=2.8, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 15, if_plot=True)

data_process.find_PSF(radius = 30, user_option = True)
data_process.plot_overview(label = 'Example', target_label = None)


#Start to produce the class and params for lens fitting.
from galight.fitting_specify import FittingSpecify

# data_process.apertures = []

fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1) #, fix_n_list= [[0,4],[1,1]])
# psf_error_map = np.ones_like(data_process.PSF_list[data_process.psf_id_for_fitting]) *0.01 # It is in the variance unit (std^2).
# fit_sepc.prepare_fitting_seq(point_source_num = 1, psf_error_map = psf_error_map)

fit_sepc.build_fitting_seq()

#Plot the initial settings for fittings. 
fit_sepc.plot_fitting_sets()

#Setting the fitting method and run.
from galight.fitting_process import FittingProcess
fit_run = FittingProcess(fit_sepc, savename = 'savename', fitting_level='norm')
fit_run.run(algorithm_list = ['PSO', 'PSO'], setting_list=[None,None])
            # setting_list = [{'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}, {'n_burn': 100, 'n_run': 100, 'walkerRatio': 10,'sigma_scale': .1}])
fit_run.plot_final_galaxy_fit()
# fit_run.dump_result()
print(fit_run.final_result_galaxy[0])
print(fit_run.final_result_ps[0])

#%%
# fig, ax = plt.subplots(figsize=(15,15))
# plt.imshow(data, norm=LogNorm(), origin='lower', cmap = my_cmap) 
# for i in range(len(pos[0])):
#     # print(pos[0][i],pos[1][i])
#     plt.scatter(pos[0][i],pos[1][i],s=30,c='k' )
# plt.show()     


# #%%

# filters = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M', 'F164N', 'F182M',
#            'F187N', 'F200W', 'F210M', 'F212N', 'F250M', 'F277W', 'F300M', 'F322W2', 'F323N',
#            'F335M', 'F356W', 'F360M', 'F405N', 'F410M', 'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']

# psf_fwhm = [0.987, 1.103, 1.298, 1.553, 1.628, 1.770, 1.801, 1.494, 1.990, 2.060, 2.141, 2.304, 2.341, 1.340,
#             1.444, 1.585, 1.547, 1.711, 1.760, 1.830, 1.901, 2.165, 2.179, 2.300, 2.302, 2.459, 2.507, 2.535, 2.574]

# zp_modA = [25.7977, 25.9686, 25.8419, 24.8878, 27.0048, 25.6536, 24.6957, 22.3073, 24.8258, 22.1775, 25.3677, 24.3296,
#            22.1036, 22.7850, 23.5964, 24.8239, 23.6452, 25.3648, 20.8604, 23.5873, 24.3778, 23.4778, 20.5588,
#            23.2749, 22.3584, 23.9731, 21.9502, 20.0428, 19.8869, 21.9002]

# zp_modB = [25.7568, 25.9771, 25.8041, 24.8738, 26.9821, 25.6279, 24.6767, 22.2903, 24.8042, 22.1499, 25.3391, 24.2909,
#            22.0574, 22.7596, 23.5011, 24.6792, 23.5769, 25.3455, 20.8631, 23.4885, 24.3883, 23.4555, 20.7007,
#            23.2763, 22.4677, 24.1562, 22.0422, 20.1430, 20.0173, 22.4086]

# dict_utils = {filters[i]: {'psf fwhm': psf_fwhm[i], 'VegaMAG zp modA': zp_modA[i],
#                            'VegaMAG zp modB': zp_modB[i]} for i in range(len(filters))}




# from photutils.detection import DAOStarFinder

# print("Using 2D Background")
# from astropy.stats import sigma_clipped_stats, SigmaClip
# sigma_clip = SigmaClip(sigma=3.)
# coverage_mask = (data == 0)
# from photutils.background import MMMBackground, MADStdBackgroundRMS, Background2D

# th = [10, 10]  # threshold level for the two filters (length must match number of filters analyzed)
# def find_stars(img, filt='F070W', var_bkg=False):
#     bkgrms = MADStdBackgroundRMS()
#     mmm_bkg = MMMBackground()

#     # image = fits.open(dict_images[det][filt]['images'][i])
#     data = img
#     print("")
#     sigma_psf = dict_utils[filt]['psf fwhm']

#     print("FWHM for the filter {f}:".format(f=filt), sigma_psf, "px")

#     if var_bkg:
#         print("Using 2D Background")
#         sigma_clip = SigmaClip(sigma=3.)
#         coverage_mask = (data == 0)

#         bkg = Background2D(data, (25, 25), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=mmm_bkg,
#                            coverage_mask=coverage_mask, fill_value=0.0)

#         data_bkgsub = data.copy()
#         data_bkgsub = data_bkgsub - bkg.background

#         _, _, std = sigma_clipped_stats(data_bkgsub)

#     else:

#         std = bkgrms(data)
#         bkg = mmm_bkg(data)

#         data_bkgsub = data.copy()
#         data_bkgsub -= bkg
#     print(std)
#     daofind = DAOStarFinder(threshold=3, fwhm=sigma_psf, roundhi=1.0,
#                             roundlo=-1.0, sharplo=0.30, sharphi=1.40)

#     found_stars = daofind(data_bkgsub)
#     # dict_aper[det][filt]['sources found'] = found_stars
#     # print("")
#     # print("Number of sources found in the image:", len(found_stars))
#     # print("-------------------------------------")
#     # print("")
#     return found_stars
# res = find_stars(data, filt=img_filter)
# from matplotlib.colors import LogNorm
# fig, ax = plt.subplots(figsize=(15,15))
# plt.imshow(data, norm=LogNorm(), origin='lower', cmap = my_cmap) 
# for j in range(len(res)):
#     plt.scatter(res[j]['xcentroid'], res[j]['ycentroid'],s=30,c='k' )
# plt.show()     
