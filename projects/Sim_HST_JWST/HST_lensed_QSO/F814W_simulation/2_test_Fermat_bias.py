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
    #file_type = 'model_result.pkl'

# folder_type = 'simulations_700_subg30/sim_lens_ID_subg30_'
# file_type = 'result_subg3_addmask.pkl'

folder_type = 'simulations_700_subg30/sim_lens_noqso_ID_subg30_'
file_type = 'result_subg3_addmask.pkl'
# file_type = 'model_result_noTD_subg3.pkl'

import glob
folder_list = glob.glob(folder_type+'*')
folder_list.sort()
test_numer = 30 #len(folder_list)
folder_list = folder_list[:test_numer]

for folder in folder_list:    
    #==============================================================================
    # ##### lens mass model 
    #==============================================================================
    from lenstronomy.LensModel.lens_model import LensModel
    lens_model_list = ['PEMD','SHEAR']
    lens_model_class = LensModel(lens_model_list)
#    print(folder)
    #Load the true parameter:
    model_lists, para_s, lens_info= pickle.load(open(folder+'/sim_kwargs.pkl','rb'))  
    kwargs_lens_list = para_s[0]
    if len(para_s[-1]['ra_image']) < 4:
#        print(para_s[-1]['ra_image'])
        para_s[-1]['ra_image'], para_s[-1]['dec_image'] = para_s[-1]['ra_image'][:2], para_s[-1]['dec_image'][:2]
    x_image, y_image = para_s[-1]['ra_image'], para_s[-1]['dec_image']
    source_pos = [para_s[2][0]['center_x'], para_s[2][0]['center_y']]
    potential = lens_model_class.potential(x_image, y_image, kwargs=kwargs_lens_list)
    fer_pot = lens_model_class.fermat_potential(x_image, y_image, kwargs_lens_list)
    geometry = potential + fer_pot
    # geometry = ((x_image - source_pos[0])**2 + (y_image - source_pos[1])**2) / 2.
#    print("potential:", potential, "geometry:", geometry)
#    
    #Load the inferred parameter:    
    multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(folder+'/'+file_type,'rb'))
    print(kwargs_result['kwargs_lens'][0])
    # x_image_infe, y_image_infe = x_image, y_image
    # kwargs_lens_list_infe = kwargs_lens_list
    
    
    kwargs_lens_list_infe = kwargs_result['kwargs_lens']
    if kwargs_result['kwargs_ps'] != []:
        x_image_infe, y_image_infe = kwargs_result['kwargs_ps'][0]['ra_image'], kwargs_result['kwargs_ps'][0]['dec_image']
    else:    
        ID = folder[-3:]
        qso_folder = folder_type[:-16] + 'ID_subg30_{0}/'.format(ID)
        _, _, kwargs_result_withQSO, _, _, _ = pickle.load(open(qso_folder+'model_result_calNoiseMap_modNoisemap_boostPossionx3_subg3.pkl','rb'))
        x_image_infe, y_image_infe = kwargs_result_withQSO['kwargs_ps'][0]['ra_image'], kwargs_result_withQSO['kwargs_ps'][0]['dec_image']
    
    source_pos_infe = [kwargs_result['kwargs_source'][0]['center_x'], kwargs_result['kwargs_source'][0]['center_y']]
    potential_infe = lens_model_class.potential(x_image_infe, y_image_infe, kwargs=kwargs_lens_list_infe)
    # geometry_infe = ((x_image_infe - source_pos_infe[0])**2 + (y_image_infe - source_pos_infe[1])**2) / 2.
    fer_pot_infe = lens_model_class.fermat_potential(x_image_infe, y_image_infe, kwargs_lens_list_infe)
    geometry_infe = potential_infe + fer_pot_infe
#    print(round(source_pos_infe[0],4), round(source_pos_infe[1],4))
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
    # potential_mismatch_list.append(np.average( true_potential_diff - infe_potential_diff ))
    # geometry_mismatch_list.append(np.average( true_geometry_diff - infe_geometry_diff ))

plt.figure(figsize=(8,6))
high0, x0, _ = plt.hist(potential_mismatch_list,normed=True, histtype=u'step',
         label=('potential mismatch'), linewidth = 2, color='orange')
high1, x1, _ = plt.hist(geometry_mismatch_list,normed=True, histtype=u'step',
         label=('geometry mismatch'), linewidth = 2, color='green')
plt.xlabel("mismatch",fontsize=27)
plt.ylabel("Density",fontsize=27)
plt.xlim(-0.02, 0.02)
plt.ylim(0,300)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.show()
