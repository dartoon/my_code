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
import astropy.io.fits as pyfits
from matplotlib.colors import LogNorm
import copy
from lenstronomy.Plots.model_plot import ModelPlot

ID_list = ['084710.40-001302.6' ,'121405.12+010205.1', '141637.44+003352.2','220906.91+004543.9',
      '233713.66+005610.8','021930.51-055643.0','022105.64-044101.5']
band_list = ['G', 'R', 'I', 'Z', 'Y']

pix_scale = 0.167

for ID in [ID_list[6]]:
    count = 0
    for band in band_list:
        picklename = '{0}/fit_image_{0}_HSC-{1}.pkl'.format(ID, band)
        result = pickle.load(open(picklename,'rb'))
        best_fit, chain_list_result, trans_paras, material = result
        source_result, image_host, ps_result, image_ps, _ =best_fit
        chain_list, _ = chain_list_result
        sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[1]
        
        if len(material) == 8:
            multi_band_list, kwargs_model, kwargs_result, QSO_msk, kwargs_fixed_source, kwargs_fixed_ps, kwargs_constraints, kwargs_numerics = material
        elif len(material) == 9:
            multi_band_list, kwargs_model, kwargs_result, QSO_msk, kwargs_fixed_source, kwargs_fixed_ps, kwargs_constraints, kwargs_numerics, classes = material

        if count == 0:
            QSO_img = multi_band_list[0][0]['image_data']
            plt.imshow(QSO_img, origin='low', norm=LogNorm())
            for i in range(len(source_result)):
                obj_x, obj_y = len(QSO_img)/2 - source_result[i]['center_x']/pix_scale, len(QSO_img)/2+source_result[i]['center_y']/pix_scale
                plt.text(obj_x, obj_y, "obj{0}".format(i), fontsize=15, color='k')
            plt.show()   
        modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result,
                                  arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[QSO_msk])
        f, axes = plt.subplots(3, 3, figsize=(16, 16), sharex=False, sharey=False)
        modelPlot.data_plot(ax=axes[0,0], text="Data")
        modelPlot.model_plot(ax=axes[0,1])
        modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
        
        modelPlot.decomposition_plot(ax=axes[1,0], text='Host galaxy', source_add=True, unconvolved=True)
        modelPlot.decomposition_plot(ax=axes[1,1], text='Host galaxy convolved', source_add=True)
        modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
        
        modelPlot.subtract_from_data_plot(ax=axes[2,0], text='Data - Point Source', point_source_add=True)
        modelPlot.subtract_from_data_plot(ax=axes[2,1], text='Data - host galaxy', source_add=True)
        modelPlot.subtract_from_data_plot(ax=axes[2,2], text='Data - host galaxy - Point Source', source_add=True, point_source_add=True)
        f.tight_layout()
        plt.show()
        count += 1