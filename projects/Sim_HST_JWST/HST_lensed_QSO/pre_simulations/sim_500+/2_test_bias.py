#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 21:20:09 2017

@author: dxh
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle

potential_mismatch_list = []
geometry_mismatch_list = []

#folder_type = 'sim_lens_ID_'
#file_type = '2nd_model_result_newlist.pkl'

folder_type = 'sim_lens_noqso_ID_'
file_type = '2nd_model_result_improve.pkl'

for seed in range(501, 522):
    #==============================================================================
    # ##### lens mass model 
    #==============================================================================
    from lenstronomy.LensModel.lens_model import LensModel
    lens_model_list = ['SPEMD','SHEAR']
    lens_model_class = LensModel(lens_model_list)
    folder = folder_type + '{0}/'.format(seed)
#    print(folder)
    #Load the true parameter:
    model_lists, para_s, lens_info= pickle.load(open(folder+'sim_kwargs.pkl','rb'))  
    kwargs_lens_list = para_s[0]
    if len(para_s[-1]['ra_image']) < 4:
#        print(para_s[-1]['ra_image'])
        para_s[-1]['ra_image'], para_s[-1]['dec_image'] = para_s[-1]['ra_image'][:2], para_s[-1]['dec_image'][:2]
    x_image, y_image = para_s[-1]['ra_image'], para_s[-1]['dec_image']
    source_pos = [para_s[2][0]['center_x'], para_s[2][0]['center_y']]
    potential = lens_model_class.potential(x_image, y_image, kwargs=kwargs_lens_list)
    geometry = ((x_image - source_pos[0])**2 + (y_image - source_pos[1])**2) / 2.
#    print("potential:", potential, "geometry:", geometry)
#    
    #Load the inferred parameter:    
    multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(folder+file_type,'rb'))
    kwargs_lens_list_infe = kwargs_result['kwargs_lens']
    x_image_infe, y_image_infe = kwargs_result['kwargs_ps'][0]['ra_image'], kwargs_result['kwargs_ps'][0]['dec_image']
    source_pos_infe = [kwargs_result['kwargs_source'][0]['center_x'], kwargs_result['kwargs_source'][0]['center_y']]
    potential_infe = lens_model_class.potential(x_image_infe, y_image_infe, kwargs=kwargs_lens_list_infe)
    geometry_infe = ((x_image_infe - source_pos_infe[0])**2 + (y_image_infe - source_pos_infe[1])**2) / 2.
    print(round(source_pos_infe[0],4), round(source_pos_infe[1],4))
#    print("potential:", potential_infe, "geometry:", geometry_inf)
    true_potential_diff = potential[0]- potential[1:]
    true_geometry_diff = geometry[0]- geometry[1:]
    infe_potential_diff = potential_infe[0]- potential_infe[1:]
    infe_geometry_diff = geometry_infe[0]- geometry_infe[1:]
#    print("A-BCD true: potential, geometry:", true_potential_diff, true_geometry_diff)
#    print("A-BCD infer: potential, geometry:",infe_potential_diff, infe_geometry_diff)
#    print("Mismatch:", true_potential_diff- infe_potential_diff, true_geometry_diff- infe_geometry_diff) 
    for i in range(len(true_potential_diff)):
        potential_mismatch_list.append(true_potential_diff[i]- infe_potential_diff[i])
        geometry_mismatch_list.append(true_geometry_diff[i]- infe_geometry_diff[i])

plt.figure(figsize=(8,6))
high0, x0, _ = plt.hist(potential_mismatch_list,normed=True, histtype=u'step',
         label=('potential mismatch'), linewidth = 2, color='orange')
high1, x1, _ = plt.hist(geometry_mismatch_list,normed=True, histtype=u'step',
         label=('geometry mismatch'), linewidth = 2, color='green')
plt.xlabel("mismatch",fontsize=27)
plt.ylabel("Density",fontsize=27)
plt.xlim(-0.01, 0.034)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.show()
