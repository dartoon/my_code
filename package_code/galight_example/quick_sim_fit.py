#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 12:32:54 2021

@author: Dartoon
"""

# some standard python imports #
import numpy as np
import matplotlib.pyplot as plt

import lenstronomy.Util.simulation_util as sim_util
import lenstronomy.Util.image_util as image_util
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.Data.psf import PSF

# data specifics
background_rms = .05  # background noise per pixel
exp_time = 100  # exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
numPix = 50  # cutout pixel size
deltaPix = 0.3  # pixel size in arcsec (area per pixel = deltaPix**2)
fwhm = 0.8  # full width half max of PSF

kwargs_data = sim_util.data_configure_simple(numPix, deltaPix, exp_time, background_rms, inverse=True)
data_class = ImageData(**kwargs_data)
# PSF specification
kwargs_psf = {'psf_type': 'GAUSSIAN', 'fwhm': fwhm, 'pixel_size': deltaPix, 'truncation': 6}
psf_class = PSF(**kwargs_psf)
# create a model with three Sersic profiles
# all the models are part of 'lens_light_model_list', meaning that their surface brightness profile are not lensed
lens_light_model_list = ['SERSIC_ELLIPSE']
from lenstronomy.LightModel.light_model import LightModel
lightModel = LightModel(lens_light_model_list)

# e1, e2 = param_util.ellipticity2phi_q(source['e1'], source['e2'])
import lenstronomy.Util.param_util as param_util
phi, q = np.pi/4, 0.5
e1, e2 = param_util.phi_q2_ellipticity(phi, q)
kwargs_1 = {'amp': 100, 'R_sersic': .5, 'n_sersic': 3, 'e1': e1, 'e2': e2, 'center_x': 0, 'center_y': 0}
kwargs_light = [kwargs_1]

# here we super-sample the resolution of some of the pixels where the surface brightness profile has a high gradient 
supersampled_indexes = np.zeros((numPix, numPix), dtype=bool)
supersampled_indexes[23:27, 23:27] = True
kwargs_numerics = {'supersampling_factor': 4, 
                   'compute_mode': 'adaptive',
                  'supersampled_indexes': supersampled_indexes}
from lenstronomy.ImSim.image_model import ImageModel
imageModel = ImageModel(data_class, psf_class, lens_light_model_class=lightModel, kwargs_numerics=kwargs_numerics)
image_sim = imageModel.image(kwargs_lens_light=kwargs_light)
poisson = image_util.add_poisson(image_sim, exp_time=exp_time)
bkg = image_util.add_background(image_sim, sigma_bkd=background_rms)
image_noisy = image_sim + bkg + poisson
plt.imshow(image_noisy, origin='lower')
plt.show()

#%%
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
data_process = DataProcess(fov_image = image_noisy, target_pos = [25, 25],
                           pos_type = 'pixel', exptime = 100,
                          rm_bkglight = False, if_plot=False, zp = 25)

data_process.PSF_list = [psf_class.kernel_pixel]

data_process.generate_target_materials(radius=25, create_mask = False, nsigma=2.8,
                                      exp_sz= 1, npixels = 10, if_plot=True, bkg_std = background_rms)
data_process.deltaPix = 1
data_process.checkout() #Check if all the materials is known.
fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 0, sersic_major_axis=False)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
fit_sepc.plot_fitting_sets()
fit_sepc.build_fitting_seq()
#Setting the fitting method and run.
fit_run = FittingProcess(fit_sepc)
fit_run.run(algorithm_list = ['PSO', 'PSO'], setting_list = [{'sigma_scale': 1., 'n_particles': 50, 'n_iterations': 50}]*2)
fit_run.plot_final_galaxy_fit()
# fit_run.dump_result()
print(phi/np.pi*180, fit_run.final_result_galaxy[0]['phi_G']/np.pi*180)
