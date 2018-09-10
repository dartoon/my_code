#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 18:04:53 2018

@author: Dartoon
"""
#import numpy as np
#import os
#import time
#import corner
#import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#
#
#from lenstronomy.SimulationAPI.simulations import Simulation
#SimAPI = Simulation()
#
#
## data specifics
#background_rms = 0.1  #  background noise per pixel (Gaussian)
#exp_time = 100.  #  exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
#numPix = 80  #  cutout pixel size
#kernel_size = 91
#deltaPix = 1  #  pixel size in arcsec (area per pixel = deltaPix**2)
#fwhm = 2.5  # full width half max of PSF
#psf_type = 'GAUSSIAN'  # 'GAUSSIAN', 'PIXEL', 'NONE'
#
#
## initial input simulation
## generate the coordinate grid and image properties
#
#data_class = SimAPI.data_configure(numPix, deltaPix, exp_time, background_rms)
## generate the psf variables
#
#psf_class = SimAPI.psf_configure(psf_type=psf_type, fwhm=fwhm, kernelsize=51,
#                                 deltaPix=deltaPix, truncate=3)
#center_x = 0.
#center_y = 0.
#
#from lenstronomy.PointSource.point_source import PointSource
#point_source_list = ['UNLENSED']
#pointSource = PointSource(point_source_type_list=point_source_list)
#kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [1]}]
#
#kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
#
#from lenstronomy.ImSim.image_model import ImageModel
#imageModel_ps = ImageModel(data_class, psf_class, point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
#image_ps = imageModel_ps.image(kwargs_ps=kwargs_ps)
#plt.matshow(np.log10(image_ps), origin='lower')
#plt.show()
#
import numpy as np
import os
import time
import copy
import corner
import astropy.io.fits as pyfits

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from lenstronomy.SimulationAPI.simulations import Simulation
SimAPI = Simulation()

# import PSF file (here as a relative path in the lenstronomy_extension repository)
# the psf_example.fits file can be found here:
# https://github.com/sibirrer/lenstronomy_extensions/tree/master/Data/PSF_TinyTim
# and imported from a local file path as well
# data specifics
background_rms = 0.1  #  background noise per pixel (Gaussian)
exp_time = 100.  #  exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
numPix = 21  #  cutout pixel size
deltaPix = 1  #  pixel size in arcsec (area per pixel = deltaPix**2)
fwhm = 2.5  # full width half max of PSF (only valid when psf_type='gaussian')
psf_type = 'GAUSSIAN'  # 'gaussian', 'pixel', 'NONE'
kernel_size = 21

# initial input simulation

# generate the coordinate grid and image properties
data_class = SimAPI.data_configure(numPix, deltaPix, exp_time, background_rms)
# generate the psf variables
psf_class = SimAPI.psf_configure(psf_type=psf_type, fwhm=fwhm, kernelsize=kernel_size, deltaPix=deltaPix, kernel=None)

# quasar center (we chose a off-centered position)
center_x = 0.07
center_y = 0.01

# quasar brightness (as measured as the sum of pixel values)
point_amp = 1
from lenstronomy.PointSource.point_source import PointSource
point_source_list = ['UNLENSED']
pointSource = PointSource(point_source_type_list=point_source_list)
kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]

from lenstronomy.LightModel.light_model import LightModel
light_model_list = ['SERSIC_ELLIPSE', 'SERSIC']
lightModel = LightModel(light_model_list=light_model_list)
import lenstronomy.Util.param_util as param_util
e1, e2 = param_util.phi_q2_ellipticity(phi=0.3, q=0.6)
kwargs_disk = {'amp': 1, 'n_sersic': 1, 'R_sersic': 0.7, 'e1': e1, 'e2': e2, 'center_x': center_x, 'center_y': center_y}
kwargs_buldge = {'amp': 1, 'n_sersic': 4, 'R_sersic': 0.3, 'center_x': center_x, 'center_y': center_y}
kwargs_host = [kwargs_disk, kwargs_buldge]

from lenstronomy.ImSim.image_model import ImageModel
kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
imageModel = ImageModel(data_class, psf_class, source_model_class=lightModel,
                                point_source_class=pointSource, kwargs_numerics=kwargs_numerics)

# simulate image with the parameters we have defined above #
image = imageModel.image(kwargs_source=kwargs_host, kwargs_ps=kwargs_ps)
imageModel_ps = ImageModel(data_class, psf_class, point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
image_ps = imageModel_ps.image(kwargs_ps=kwargs_ps)
plt.matshow((image_ps), origin='lower')
plt.show()