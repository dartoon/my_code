#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 14:45:16 2020

@author: Xuheng Ding
"""

#photutils in version > = 0.7.2
#astropy in version astropy-4.0.1

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

#%%
fitsFile = pyfits.open('../example_data/HSC/QSO/000017.88+002612.6_HSC-I.fits')

fov_image= fitsFile[1].data
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
err_data= fitsFile[3].data ** 0.5

file_header0 = fitsFile[0].header
FLUXMAG0 = file_header0['FLUXMAG0']
zp =  2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
PSF = pyfits.getdata('../example_data/HSC/QSO/000017.88+002612.6_HSC-I_psf.fits')

#%%Start to use galight
from galight.data_process import DataProcess
QSO_RA = 0.07452999800443649
QSO_DEC = 0.4368380010128021
data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [QSO_RA, QSO_DEC],
                           pos_type = 'wcs', header = header,
                          rm_bkglight = True, if_plot=False, zp = zp)

data_process.generate_target_materials(radius=None, create_mask = False, nsigma=2.8,
                                      exp_sz= 1.2, npixels = 15, if_plot=True, use_moments=False)
data_process.PSF_list = [PSF]

data_process.checkout() #Check if all the materials is known.

# %%
from galight.tools.measure_tools import image_moments
image = data_process.target_stamp
segm = data_process.segm_deblend
moments = image_moments(image, segm, 1)

target = data_process.apertures[0]
positions = [moments['X'],moments['Y']]
a = moments['Mrr']
b = moments['Mrr']*moments['q']
theta = moments['phi_deg']*np.pi/180



#%%Update the apretures:
import copy
apertures = copy.deepcopy(data_process.apertures)
add_aperture0 = copy.deepcopy(apertures[0])
add_aperture0.a, add_aperture0.b = 3, 3 # A a b for aperture
add_aperture0.positions = np.array([40,50])
data_process.apertures = apertures +  [add_aperture0]#Pass apertures to the data_process

# #Update the segm map:
from galight.tools.measure_tools import mask_obj
mask = mask_obj(data_process.target_stamp, [add_aperture0])
data_process.segm_deblend = data_process.segm_deblend*mask[0] + (1-mask[0])*(np.max(data_process.segm_deblend)+1)
plt.imshow(data_process.segm_deblend, origin='lower')
plt.colorbar()
plt.show()

# %%Start to produce the class and params for lens fitting.
from galight.fitting_specify import FittingSpecify
fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor=3, mpi=False)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
# fit_sepc.plot_fitting_sets()
fit_sepc.build_fitting_seq()

#Setting the fitting method and run.
from galight.fitting_process import FittingProcess
fit_run = FittingProcess(fit_sepc, savename = 'HSC_QSO', fitting_level='norm')
fit_run.run(algorithm_list = ['PSO','PSO'], setting_list=None,threadCount=2)
# fit_run.plot_all()
# fit_run.dump_result()
print(fit_run.final_result_galaxy[0])

# #%%
# import pickle
# from galight.tools.asymmetry_tools import Measure_asy
# from galight.tools.astro_tools import plt_fits

# fit_run_pkl = fit_run
# obj_id = 0
# asy_class = Measure_asy(fit_run_pkl, seg_cal_reg = 'or', obj_id=obj_id, #consider_petrosian=True, 
#                         extend=1.5, eta=0.2)

# plt_fits(asy_class.img,colorbar=True)

# asy_class.asy_segm(mask_type='segm')
# # asy_class.asy_segm(mask_type='aper')  #!!!
# pos = asy_class.find_pos()
# print(pos["x"])
# asy_class.run_bkg(rotate_pix = pos["x"], if_remeasure_bkg=True)
# asy = asy_class.cal_asymmetry(rotate_pix = pos["x"],  #!!!
#                               if_plot=True, if_plot_bkg=True)

# from galight.tools.asymmetry_tools import cal_r_petrosian
# r_p = cal_r_petrosian(asy_class.img,
#                 center = np.array([len(asy_class.img)/2]*2)+pos['x'],
#                 if_plot=True)

# # plt.imshow(asy_class._segm, origin='lower')
# # plt.show()
# #%%Test for write cal_r_petrosian
# from galight.tools.measure_tools import SB_profile
# # image = fit_run.fitting_specify_class.kwargs_data['image_data']
# # mask = (asy_class.segm== asy_class.segm_id)
# # apertures = data_process.apertures
# # q, theta = apertures[obj_id].b/apertures[obj_id].a, apertures[obj_id].theta
# # r_SB, r_grids  =  SB_profile(data_process.target_stamp*mask, 
# #                               center = asy_class.apertures[asy_class.obj_id].positions, 
# #                               if_plot=True, fits_plot = True, if_annuli= False, grids=40, q =q, theta=theta)
# # print(r_grids.max())
# # r_SB_annu, _  =  SB_profile(data_process.target_stamp*mask, 
# #                             center = asy_class.apertures[asy_class.obj_id].positions, 
# #                               if_plot=False, fits_plot = True, if_annuli= True, grids=40,  q =q, theta=theta)
# # r_p = r_grids[np.sum(r_SB_annu/r_SB>0.2)]
# # print(r_p)
# # from galight.tools.measure_tools import mask_obj
# # aper = asy_class.apertures[asy_class.obj_id]
# # aper.a, aper.b = r_p*1.5, r_p*1.5
# # mask = mask_obj(data_process.target_stamp, [aper])[0]
# # # plt.imshow(mask)

# r_p = cal_r_petrosian(asy_class.img,
#                       center = np.array([len(asy_class.img)/2]*2)+pos['x'],
#                         q = data_process.apertures[obj_id].b/data_process.apertures[obj_id].a,
#                         theta = data_process.apertures[obj_id].theta,
#                         # q = 1,
#                         # theta = 0,
#                       if_plot=True)

# from galight.tools.cutout_tools import pix_region
# region = pix_region(pos['x']+np.array([len(asy_class.img)/2]*2),
#                     r_p, q=data_process.apertures[obj_id].b/data_process.apertures[obj_id].a,
#                       theta = data_process.apertures[obj_id].theta)
# mask = region.to_mask(mode='exact')

# from photutils import EllipticalAperture
# ap = EllipticalAperture( pos['x']+np.array([len(asy_class.img)/2]*2), 
#                         a = r_p, 
#                         b = r_p * data_process.apertures[obj_id].b/data_process.apertures[obj_id].a,
#                         theta=data_process.apertures[obj_id].theta)

# from galight.tools.measure_tools import mask_obj
# mask_ap = mask_obj(asy_class.img, [ap])[0]

# data = mask.cutout(asy_class.img)
# plt.imshow(mask.data*data)
# print(np.sum(mask.data*data))
# print(np.sum((1-mask_ap)*asy_class.img))
# plt.imshow(mask_ap, origin='lower')
# plt.show()
# # print(r_p)
