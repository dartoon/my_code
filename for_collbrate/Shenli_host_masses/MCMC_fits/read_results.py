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

ID_list = ['084710.40-001302.6' ,'121405.12+010205.1', '141637.44+003352.2','220906.91+004543.9',
      '233713.66+005610.8','021930.51-055643.0','022105.64-044101.5']
band_list = ['G', 'R', 'I', 'Z', 'Y']
which_ID = input("To check which ID?\n0: {0}\n1: {1}'\n2: {2}\n3: {3}\n4: {4}\n5: {5}\n6: {6}\n".format(ID_list[0], ID_list[1], ID_list[2], ID_list[3], ID_list[4], ID_list[5],ID_list[6]))
which_band = input("select a band: \n0:'G' \n1:'R'\n2:'I'\n3:'Z'\n4:'Y'\n")
which_ID, which_band = int(which_ID), int(which_band)

ID = ['084710.40-001302.6' ,'121405.12+010205.1', '141637.44+003352.2','220906.91+004543.9',
      '233713.66+005610.8','021930.51-055643.0','022105.64-044101.5'][which_ID]
band = ['G', 'R', 'I', 'Z', 'Y'][which_band]
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

# print("best-fit source_result:", source_result)
# print("best-fit ps_result:", ps_result)

#If want to see the corner map for all the fitting parameters:
# plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
# plt.show()


pix_scale = 0.167
QSO_img = multi_band_list[0][0]['image_data']    

#If need to translate:
# from lenstronomy.Sampling.parameters import Param
# from lenstronomy.ImSim.image_linear_solve import ImageLinearFit
# param = Param(kwargs_model, kwargs_fixed_source=kwargs_fixed_source, **kwargs_constraints)
# kwargs_out = param.args2kwargs(samples_mcmc[0])
# kwargs_light_source_out = kwargs_out['kwargs_source']
# kwargs_ps_out =  kwargs_out['kwargs_ps']
# data_class, psf_class, lightModel, pointSource = classes
# imageLinearFit = ImageLinearFit(data_class=data_class, psf_class=psf_class,
#                                 source_model_class=lightModel,
#                                 point_source_class=pointSource, 
#                                 kwargs_numerics=kwargs_numerics)   
# image_reconstructed, _, _, _ = imageLinearFit.image_linear_solve(kwargs_source=kwargs_light_source_out, kwargs_ps=kwargs_ps_out)
# from lenstronomy.ImSim.image_model import ImageModel
# imageModel = ImageModel(data_class, psf_class, source_model_class=lightModel,
#                         point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
# image_j = imageModel.source_surface_brightness(kwargs_light_source_out,unconvolved= False, k=1)    

#%% Recover the plot
from lenstronomy.Plots.model_plot import ModelPlot
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

import sys
sys.path.insert(0, '../../../py_tools/')
from transfer_to_result import transfer_to_result
if len(image_ps) >1:
    image_ps = np.sum(image_ps, axis=0)
else:
    image_ps = image_ps[0]
_ = transfer_to_result(data=QSO_img, pix_sz = pix_scale,  
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, error_map=multi_band_list[0][0]['noise_map'],
                            zp=27.0, fixcenter=True,ID='Data', QSO_msk = QSO_msk, tag=None, plot_compare = True)


#%%Recover the translated cornor plot
plt.imshow(QSO_img, origin='low', norm=LogNorm())
for i in range(len(ps_result)):
    obj_x, obj_y = len(QSO_img)/2 - ps_result[i]['ra_image'][0]/pix_scale, len(QSO_img)/2+ps_result[i]['dec_image'][0]/pix_scale
    # print(obj_x, obj_y)
    plt.text(obj_x, obj_y, "QSO{0}".format(i), fontsize=10, color='k')
plt.show()    

plt.imshow(QSO_img, origin='low', norm=LogNorm())
for i in range(len(source_result)):
    obj_x, obj_y = len(QSO_img)/2 - source_result[i]['center_x']/pix_scale, len(QSO_img)/2+source_result[i]['center_y']/pix_scale
    plt.text(obj_x, obj_y, "obj{0}".format(i), fontsize=15, color='k')
plt.show()   


#%%The host flux for the host flux:
#Read the fitting parameter.
idx = 3  # this values can be changed to check the inference for difference parameter.
# #Plot the translated host MCMC corner inference for the QSO and host inference.  
v_l=np.percentile(samples_mcmc[:,idx],16,axis=0)
v_m=np.percentile(samples_mcmc[:,idx],50,axis=0)
v_h=np.percentile(samples_mcmc[:,idx],84,axis=0)
print(param_mcmc[idx], ":", v_l, v_m, v_h)


mcmc_new_list, labels_new, _ = trans_paras
zp = 27.0
#For the translated totol flux.
mcmc_new_list = np.asarray(mcmc_new_list)
idx = 0
v_l=np.percentile(mcmc_new_list[:,idx],16,axis=0)
v_m=np.percentile(mcmc_new_list[:,idx],50,axis=0)
v_h=np.percentile(mcmc_new_list[:,idx],84,axis=0)
print(labels_new[idx], "(lower mid upper): {0:.3f} {1:.3f} {2:.3f}".format(v_l, v_m, v_h) )
#Print the magnitude inference for a given host:
mag_l, mag_m, mag_h = -2.5*np.log10(v_h) + zp,  -2.5*np.log10(v_m) + zp,  -2.5*np.log10(v_l) + zp    
print(labels_new[idx].split(' ')[0], " magnitude (lower mid upper): {0:.3f} {1:.3f} {2:.3f}".format(mag_l, mag_m, mag_h) )

#%%
ifplot = input("Do you want to plot the PSO and MCMC inference?\n1:Yes, else:No\n")
if int(ifplot) == 1:   
    #If want to see the corner map for all the fitting parameters:
    # plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
    # plt.show()
    
    plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
    plt.show()
    #diagnose the PSO chain convergency in the fitting
    # from lenstronomy.Plots import chain_plot
    # f, axes = chain_plot.plot_chain_list(chain_list,0)
    # plt.show()
    # #test the MCMC chain convergency in the fitting
    import lenstronomy.Plots.chain_plot as out_plot
    fig = plt.figure(figsize=(20, 15))
    ax = fig.add_subplot(111)
    out_plot.plot_mcmc_behaviour(ax, samples_mcmc, param_mcmc, dist_mcmc)   
    plt.show()        
