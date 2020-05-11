#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 21:20:09 2017

@author: dxh
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import astropy.io.fits as pyfits
import copy
import time
from lenstronomy.Util import constants as const
import lenstronomy.Util.param_util as param_util
import glob
#file name:
filt='f160w'
import pickle
import sys
sys.path.insert(0,'../../../py_tools/')
from flux_profile import cr_mask
from lenstronomy.Plots.model_plot import ModelPlot
from mask_objects import find_loc_max


folder_type = 'simulations_700_sim_subg3/sim_lens_noqso_ID_'
# # file_type = 'model_result_noTD.pkl'
file_type_list = ['model_result_noTD.pkl', 'model_result_noTD_calNoisemap_subg3.pkl', 'model_result_noTD_calNoisemap_subg4.pkl', 'model_result_noTD_subg5.pkl']
# file_type_list = ['model_result_noTD_subg3.pkl', 'model_result_noTD_calNoisemap_subg3.pkl','model_result_noTD_calNoisemap_subg3_boostnoisemap.pkl']

# folder_type = 'sim_lens_noqso_nodrz_ID_'
# file_type_list = ['model_result_noTD_subg3.pkl']

folder_list = glob.glob(folder_type+'*')
folder_list.sort()
test_numer = 30
folder_list = folder_list[:test_numer]

# del folder_list[14]
# del folder_list[13]
# del folder_list[12]

id_range = int(folder_list[0][-3:]), int(folder_list[-1][-3:])

result_dic = {}
# if 'noTD' not in file_type:
for folder in folder_list:
    folder = folder+'/'
    for file_type in file_type_list:
        read_file = folder+ file_type
        print(read_file)
        model_lists, para_s, lens_info= pickle.load(open(folder+'sim_kwargs.pkl','rb'))
        lens_model_list, lens_light_model_list, source_model_list, point_source_list = model_lists
        z_l, z_s, TD_distance, TD_true, TD_obs, TD_err_l = lens_info
        kwargs_lens_list, kwargs_lens_light_list, kwargs_source_list, kwargs_ps = para_s
        solver_type = 'PROFILE_SHEAR'
        if len(kwargs_ps['ra_image']) <4:
            kwargs_ps['ra_image'] = kwargs_ps['ra_image'][:2] 
            kwargs_ps['dec_image'] = kwargs_ps['dec_image'][:2]
            kwargs_ps['point_amp'] = kwargs_ps['point_amp'][:2]
            TD_obs = TD_obs[:2]
            TD_err_l = TD_err_l[:2]
            solver_type = 'THETA_E_PHI'
        kwargs_constraints = {'joint_source_with_point_source': [[0, 0]],
                              'num_point_source_list': [len(kwargs_ps['ra_image'])],
                              'solver_type': solver_type,  # 'PROFILE', 'PROFILE_SHEAR', 'ELLIPSE', 'CENTER'
                              'Ddt_sampling': True,
                                      }
    #    if glob.glob(folder+'model_result.pkl') == []:
    #        result_dic[folder] = [None, None]
    #        continue
        multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(read_file,'rb'))
        fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
        
        modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat")
        logL = modelPlot._imageModel.likelihood_data_given_model(source_marg=False, linear_prior=None, **kwargs_result)
        n_data = modelPlot._imageModel.num_data_evaluate
        chisq = -logL * 2 / n_data
        
        truth_dic = {}
        truth_dic['kwargs_lens'] =kwargs_lens_list
        truth_dic['kwargs_source'] =kwargs_source_list
        truth_dic['kwargs_lens_light'] =kwargs_lens_light_list
        truth_dic['kwargs_ps'] = kwargs_ps
        truth_dic['D_dt'] = TD_distance
        result_dic[folder[:-1]+file_type] = [truth_dic, kwargs_result, chisq]
        
    
#%%test parameter bias:
#para = 'theta_E'  #'gamma'
#para = 'gamma'
which = ['kwargs_lens', 'gamma']
#which = ['kwargs_source', 'center_y']
fig, ax = plt.subplots(figsize=(11,8))
c_list = ['red', 'blue', 'green', 'black']
gamma_bias_list = []
for folder in folder_list:
    ID = folder[-3:]
    read_list = [1]
    for i in read_list:
        file_type = file_type_list[i]
        key = folder_type + '{0}'.format(ID) + file_type
        gamma_bias = result_dic[key][1][which[0]][0][which[1]] - result_dic[key][0][which[0]][0][which[1]]
        gamma_bias_list.append(gamma_bias)
        plt.plot(np.linspace(0, 37), np.linspace(0, 37)*0)
        if file_type_list[i]  == 'model_result_noTD.pkl':
            label = 'model_result_noTD_subg2.pkl'
        else:
            label = file_type_list[i]
        plt.scatter(int(ID)-700, gamma_bias,
                    c=c_list[i],s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7, label = label)
        plt.text(int(ID)-700, gamma_bias, repr(round(result_dic[key][-1], 1)),fontsize=15)
#    ax.set_xticks(range(id_range[0]-1, id_range[1]+1,3)) 
#    plt.plot(np.linspace(id_range[0]-1, id_range[1]+1), np.linspace(id_range[0]-1, id_range[1]+1)*0)
    plt.xlabel("ID",fontsize=27)
    plt.ylabel(which[1]+" bias (inferred - truth)",fontsize=27)
    plt.ylim(-0.4,0.4)
    plt.tick_params(labelsize=20)
    (lines, labels) = plt.gca().get_legend_handles_labels()
    # plt.legend()
    plt.legend(lines[:len(read_list)], labels[:len(read_list)], prop={'size': 15.5, 'family': 'Arial'})
plt.show()
print(np.mean(gamma_bias_list), np.std(gamma_bias_list))


