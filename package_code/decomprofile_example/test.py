#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 16:01:17 2020

@author: Xuheng Ding

You can skip this step if the QSO stamp, noise level and the PSF is ready.
"""
#photutils in version > = 0.7.2
#astropy in version astropy-4.0.1

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

#%%
fitsFile = pyfits.open('../example_data/HST/QSO/1147_final_drz.fits')
fov_image = fitsFile[1].data # check the back grounp
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']

#%%
wht = fitsFile[2].data # The WHT map
import galight.tools.astro_tools as astro_tools
exp =  astro_tools.read_fits_exp(fitsFile[0].header)  #Read the exposure time 
mean_wht = exp * (0.0642/0.135)**2
exp_map = exp * wht/mean_wht

#%%Start to use galight
from galight.data_process import DataProcess
data_process = DataProcess(fov_image = fov_image, target_pos = [1142., 637.], pos_type = 'pixel', header = header,
                          rm_bkglight = False, exptime = exp_map, if_plot=False, zp = 27.0)
data_process.generate_target_materials(radius=65, create_mask = True, 
                                       detect_tool = 'sep' , 
                                       if_select_obj= True,
                                       nsigma=2.5,
                                       thresh = 2.5,
                                       exp_sz= 1.2, npixels = 15, if_plot=True)

#%%
data = data_process.target_stamp
from matplotlib.colors import LogNorm
from astropy.visualization import AsinhStretch, ImageNormalize, MinMaxInterval, SqrtStretch

norm = ImageNormalize(data, interval=MinMaxInterval(),
                      stretch=AsinhStretch())

plt.imshow(data, origin='lower', norm=norm, vmax = data.max(), vmin =  data.min())
plt.colorbar()
plt.show()

#%%
# from astropy.stats import gaussian_fwhm_to_sigma
# from astropy.convolution import Gaussian2DKernel
# sigma = 3.0 * gaussian_fwhm_to_sigma # FWHM = 3.
# kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
# kernel.normalize()

# import sep
# data = data_process.target_stamp
# data = data.byteswap().newbyteorder()
# objects, segment = sep.extract(data, 2.5, err=data_process.bkg_std,mask=None, minarea=5,
#                 filter_kernel=None, filter_type='matched',
#                 deblend_nthresh=32, deblend_cont=0.005, clean=True,
#                 clean_param=1.0, segmentation_map=True)


# from matplotlib.colors import LogNorm
# fig, ax = plt.subplots()
# plt.imshow(data, origin='lower', norm=LogNorm(), vmax = data.max(), vmin = 1.e-4)
# plt.colorbar()
# from matplotlib.patches import Ellipse
# # plot an ellipse for each object
# for i in range(len(objects)):
#     e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
#                 width=6*objects['a'][i],
#                 height=6*objects['b'][i],
#                 angle=objects['theta'][i] * 180. / np.pi)
#     e.set_facecolor('none')
#     e.set_edgecolor('red')
#     ax.add_artist(e)
# plt.show()  


#%%PSF works.
data_process.find_PSF(radius = 30, user_option = False)
# data_process.find_PSF(radius = 50, PSF_pos_list = [[ 350., 1055.], [2078., 1910.]], user_option = False)
# data_process.plot_overview(label = 'Example', target_label = None)
# data_process.profiles_compare(norm_pix = 5, if_annuli=False, y_log = False,
#                   prf_name_list = (['target'] + ['PSF{0}'.format(i) for i in range(len(data_process.PSF_list))]) )

# data_process.psf_id_for_fitting = 11

data_process.checkout() #Check if all the materials is known.

#%%Start to produce the class and params for lens fitting.
from galight.fitting_specify import FittingSpecify

# data_process.apertures = []

fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1) #, fix_n_list= [[0,4],[1,1]])
# psf_error_map = np.ones_like(data_process.PSF_list[data_process.psf_id_for_fitting]) *0.01 # It is in the variance unit (std^2).
# fit_sepc.prepare_fitting_seq(point_source_num = 1, psf_error_map = psf_error_map)

fit_sepc.build_fitting_seq()

#Plot the initial settings for fittings. 
fit_sepc.plot_fitting_sets()

#%%Setting the fitting method and run.
from galight.fitting_process import FittingProcess
fit_run = FittingProcess(fit_sepc, savename = 'savename', fitting_level='norm')
fit_run.run(algorithm_list = ['PSO', 'MCMC'], setting_list=[None,{'n_burn': 100, 'n_run': 100, 'walkerRatio': 10,'sigma_scale': .1}])
            # setting_list = [{'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}, {'n_burn': 100, 'n_run': 100, 'walkerRatio': 10,'sigma_scale': .1}])
# fit_run.plot_all()

# fit_run.dump_result()
# print(fit_run.final_result_galaxy[0])

#%%
# # Test load pkl
# import pickle
# picklename = 'savename.pkl'
# fitting_run_class = pickle.load(open(picklename,'rb'))
# fitting_run_class.plot_final_qso_fit()
# from galight.tools.plot_tools import profile_plots

# data = fitting_run_class.fitting_specify_class.kwargs_data['image_data']

# ps_list = fitting_run_class.image_ps_list
# galaxy_list = fitting_run_class.image_host_list
# ps_image = np.zeros_like(ps_list[0])
# for i in range(len(ps_list)):
#     ps_image = ps_image+ps_list[i]
# galaxy_image = np.zeros_like(data)
# for i in range(len(galaxy_list)):
#     galaxy_image = galaxy_image+galaxy_list[i]
# model = ps_image + galaxy_image

# flux_list_2d = [data, model, (data-model)/data]
# label_list_2d = ['data', 'model', '(data-model)/model']
# flux_list_1d = [data, model]
# label_list_1d = ['data', 'model']      
# profile_plots(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d,
#               deltaPix = fitting_run_class.fitting_specify_class.deltaPix,
#               target_ID =  'ID', if_annuli=True)

