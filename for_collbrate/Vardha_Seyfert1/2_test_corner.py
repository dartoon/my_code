#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 21:48:10 2019

@author: Dartoon
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner

picklename = 'dump_fitting.pkl' 
result = pickle.load(open(picklename,'rb'))
[source_result, image_host, ps_result, image_ps, samples_mcmc, param_mcmc, paras] = result

#%%
# here the (non-converged) MCMC chain of the non-linear parameters
print 'The original plot:'
[source_params_2, ps_param_2, mcmc_new_list, labels_new] = paras
plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
plt.show()

import astropy.io.fits as pyfits
psf = pyfits.getdata('PSF0.fits')  # Input the PSF0 to fit.

ID = 'l106'
agn_image = pyfits.getdata('{0}_cutout.fits'.format(ID))
fit_frame_size = 81
ct = (len(agn_image)-fit_frame_size)/2     # If want to cut to 81, agn_image[ct:-ct,ct:-ct]
agn_image = agn_image[ct:-ct,ct:-ct]

kwargs_constraints = {'num_point_source_list': [1]}

point_source_list = ['UNLENSED']
light_model_list = ['SERSIC_ELLIPSE'] * 3
kwargs_model = { 'source_light_model_list': light_model_list,
                'point_source_model_list': point_source_list
                }

from lenstronomy.Sampling.parameters import Param
param = Param(kwargs_model, kwargs_fixed_source=source_params_2, kwargs_fixed_ps=ps_param_2, **kwargs_constraints)
mcmc_new_list = []

from lenstronomy.ImSim.image_model import ImageModel

import lenstronomy.Util.simulation_util as sim_util
from lenstronomy.Data.imaging_data import Data
kwargs_data = sim_util.data_configure_simple(81, 0.04, 1000, 0.001, inverse=True)
data_class = Data(kwargs_data)
data_class.update_data(agn_image)
kernel_size = len(psf)
kernel = psf
kwargs_psf =  sim_util.psf_configure_simple(psf_type='PIXEL', fwhm=0.2, kernelsize=kernel_size, deltaPix=0.04, kernel=kernel)
from lenstronomy.Data.psf import PSF
psf_class = PSF(kwargs_psf)

from lenstronomy.LightModel.light_model import LightModel
light_model_list = ['SERSIC_ELLIPSE'] * len(source_params_2)
lightModel = LightModel(light_model_list=light_model_list)

from lenstronomy.PointSource.point_source import PointSource
pointSource = PointSource(point_source_type_list=point_source_list)
kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
#kwargs_numerics['mask'] = QSO_msk
imageModel = ImageModel(data_class, psf_class, source_model_class=lightModel,
                        point_source_class=pointSource, kwargs_numerics=kwargs_numerics)

#labels_new = ["Quasar flux"] +  ["host{0} flux".format(i) for i in range(3)] + ["host{0} Reff".format(i) for i in range(3)]
for i in range(len(samples_mcmc)/10):
    kwargs_lens_out, kwargs_light_source_out, kwargs_light_lens_out, kwargs_ps_out, kwargs_cosmo = param.args2kwargs(samples_mcmc[i+ len(samples_mcmc)/10*9])
    image_reconstructed, _, _, _ = imageModel.image_linear_solve(kwargs_source=kwargs_light_source_out, kwargs_ps=kwargs_ps_out)
    image_ps = imageModel.point_source(kwargs_ps_out)
    flux_quasar = np.sum(image_ps)
    fluxs, reffs = [],[]
    for j in range(3):
        image_j = imageModel.source_surface_brightness(kwargs_light_source_out,unconvolved= False, k=j)
        fluxs.append(np.sum(image_j))
        reffs.append(kwargs_light_source_out[j]['R_sersic'])
    mcmc_new_list.append([flux_quasar] + fluxs + reffs)
    if i/1000 > (i-1)/1000 :
        print "finished translate:", i      

print 'The new translated plot:'
plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
plt.show()