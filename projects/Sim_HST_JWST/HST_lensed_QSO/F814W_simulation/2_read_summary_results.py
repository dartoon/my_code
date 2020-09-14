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
filt='f814w'
import pickle
import sys
sys.path.insert(0,'../../../../py_tools/')
from flux_profile import cr_mask
from lenstronomy.Plots.model_plot import ModelPlot
from mask_objects import find_loc_max

result_dic = {}

# folder_list = use_folder #['simulations_700_subg30/sim_lens_noqso_ID_subg30_' + use_folder[i][-3:] for i in range(len(use_folder))]
# folder_list = ['simulations_700_subg30/sim_lens_noqso_ID_subg30_' + use_folder[i][-3:] for i in range(len(use_folder))]

# # folder_type  = folder_list[0][:-3]
# # file_type = 'model_result_use_drz_Noisemap_subg3.pkl'
# file_type = 'model_result_calNoiseMap_modNoisemap_boostPossionx3_subg3.pkl'
# # file_type = 'model_result_calNoiseMap_modNoisemap_boostPossionx3_subg3.pkl'
file_type = 'result_subg3_addmask.pkl'
folder_type = 'simulations_700_subg30/sim_lens_ID_subg30_'


# # # # file_type = 'model_result_subg3.pkl'
# # file_type = 'model_result_calNoiseMap_modNoisemap_boostPossionx8_noPSFerr_subg2_fixgamma.pkl'
# # # file_type = 'model_result_calNoiseMap_modNoisemap_useGrad_noPSFerr_subg3.pkl'
# file_type = 'result_subg3.pkl'
# folder_type = 'simulations_700_subg30/sim_lens_noqso_ID_subg30_'
# file_type = 'result_subg3_addmask.pkl'
# folder_type = 'simulations_700_subg30/sim_lens_noqso_ID_subg30_'



folder_list = glob.glob(folder_type+'*')
folder_list.sort()
test_numer = 30 #len(folder_list)
folder_list = folder_list[:test_numer]


id_range = int(folder_list[0][-3:]), int(folder_list[-1][-3:])


for folder in folder_list:
    folder = folder+'/'
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
    mcmc_new_list = np.array(mcmc_new_list)
    H0_list = mcmc_new_list[:,-1]
    
    lens_data = pyfits.getdata(folder+'Drz_QSO_image.fits')
    lens_mask = cr_mask(lens_data, 'normal_mask.reg')
    framesize = len(multi_band_list[0][0]['image_data'])  #81
    ct = int((len(lens_data) - framesize)/2)
    lens_data = lens_data[ct:-ct,ct:-ct]
    lens_mask = (1-lens_mask)[ct:-ct,ct:-ct]
    
    x, y =find_loc_max(lens_data)
    x_s, y_s = [], []
    center = (framesize-1)/2
    deltaPix = 0.08
    for i in range(len(x)):
        x0, y0 =  (float(x[i]) - center) * deltaPix , (float(y[i]) - center) * deltaPix
#            print(x0, y0)
        ds = (x0- kwargs_ps['ra_image'] )**2 + (y0-kwargs_ps['dec_image'])**2
        if ds.min()<0.01:
            x_s.append(x[i])
            y_s.append(y[i])
    y_grid,x_grid = np.indices((framesize,framesize))   #with higher resolution 60*6
    for i in range(len(x_s)):
        lens_mask[np.sqrt((y_grid-y_s[i])**2 + (x_grid-x_s[i])**2) <4] = 0
        
    # plt.imshow(multi_band_list[0][0]['noise_map'],origin='lower')
    # plt.show()
    # plt.imshow(lens_mask,origin='lower')
    # plt.show()        
    
    modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[lens_mask])
    logL = modelPlot._imageModel.likelihood_data_given_model(source_marg=False, linear_prior=None, **kwargs_result)
    n_data = modelPlot._imageModel.num_data_evaluate
    chisq = -logL * 2 / n_data
    if chisq< 0.0:
        print(folder[-4:-1], round(np.median(H0_list), 3))
        f, axes = modelPlot.plot_main()
        f.show()
        plt.show()
    truth_dic = {}
    truth_dic['kwargs_lens'] =kwargs_lens_list
    truth_dic['kwargs_source'] =kwargs_source_list
    truth_dic['kwargs_lens_light'] =kwargs_lens_light_list
    truth_dic['kwargs_ps'] = kwargs_ps
    truth_dic['D_dt'] = TD_distance
    result_dic[folder[:-1]] = [truth_dic, kwargs_result, 
                               [np.percentile(H0_list,16), np.percentile(H0_list,50), np.percentile(H0_list,84)], 
                               (np.percentile(chain_list[-1][1][:,0],84) - np.percentile(chain_list[-1][1][:,0],16))/2, chisq]

#%%
H0_true = 73.907
fig, ax = plt.subplots(figsize=(11,8))
H0_list = []
use_folder = []
for folder in folder_list:
    ID = folder[-3:]
    key = folder_type + '{0}'.format(ID)
    # if abs(result_dic[key][-1]) < 3.0 and result_dic[key][2][1] < 90 and result_dic[key][2][1] >61:      #Use folders meets this requirments
    if abs(result_dic[key][-1]) < 30 and result_dic[key][2][1] < 100 and result_dic[key][2][1] >60:      #Use folders meets this requirments
        ID = int(ID)
        H0 = result_dic[key][2]
        plt.scatter(ID, H0[1],
                    c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
        plt.errorbar(ID, H0[1], yerr = [[H0[2]-H0[1]], [H0[1]-H0[0]]],
                    ecolor='black', fmt='o', zorder=-500,markersize=1)  
        plt.text(ID, H0[1], repr(round(result_dic[key][-1], 1)),fontsize=15)
        H0_list.append([H0[1], (H0[2]-H0[1]+ H0[1]-H0[0])/2 ])
        use_folder.append(folder)
#        ax.set_xticks(range(id_range[0]-1, id_range[1]+1,3)) 
plt.plot(np.linspace(id_range[0]-1, id_range[1]+1), np.linspace(id_range[0]-1, id_range[1]+1)*0 + H0_true)
plt.xlabel("ID",fontsize=27)
plt.ylabel("$H_0$",fontsize=27)
plt.ylim(50,100)
plt.tick_params(labelsize=20)
plt.show()

H0_list = np.array(H0_list)
plt.hist(H0_list[:, 0])
plt.show()


submit_sum = np.array(H0_list)
print(np.mean(submit_sum[:,0]), np.std(submit_sum[:,0]))

"""
goodness_sum = round(1/float(len(submit_sum))*np.sum(((submit_sum[:,0]-H0_true)/submit_sum[:,1])**2),3)
precision_sum = round(1/float(len(submit_sum))*np.sum(submit_sum[:,1]/H0_true)*100, 3)
accuracy_sum = round(1/float(len(submit_sum))*np.sum((submit_sum[:,0]-H0_true)/H0_true)*100, 3)

# print(goodness_sum, precision_sum, accuracy_sum)

#%%test H0 Bias as function of other parameter's bias
#para = 'theta_E'  #'gamma'
#para = 'gamma'
which = ['kwargs_lens', 'gamma']
#which = ['kwargs_source', 'center_y']
fig, ax = plt.subplots(figsize=(11,8))
for folder in folder_list:
    ID = folder[-3:]
    key = folder_type + '{0}'.format(ID)
    H0 = result_dic[key][2]
    gamma_bias = result_dic[key][1][which[0]][0][which[1]] - result_dic[key][0][which[0]][0][which[1]]
    plt.scatter(gamma_bias, H0[1] - H0_true,
                c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
    plt.errorbar(gamma_bias, H0[1] - H0_true,
                 yerr = [[H0[2]-H0[1]], [H0[1]-H0[0]]],
                ecolor='black', fmt='o', zorder=-500,markersize=1)      
#    ax.set_xticks(range(id_range[0]-1, id_range[1]+1,3)) 
plt.xlabel(which[1]+" bias (inferred - truth)", fontsize=27)
plt.ylabel("$H_0$ bias (inferred - truth)", fontsize=27)
plt.ylim(-20,20)
plt.tick_params(labelsize=20)
plt.show()
"""

#%%test parameter bias:
#para = 'theta_E'  #'gamma'
#para = 'gamma'
which = ['kwargs_lens', 'gamma']
#which = ['kwargs_source', 'center_y']
fig, ax = plt.subplots(figsize=(11,8))
ct = 0
gamma_bias_list = []
for folder in folder_list:
    ID = folder[-3:]
    key = folder_type + '{0}'.format(ID)
    gamma_bias = result_dic[key][1][which[0]][0][which[1]] - result_dic[key][0][which[0]][0][which[1]]
    gamma_bias_list.append(gamma_bias)
    plt.scatter(ct, gamma_bias,
                c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
    plt.errorbar(ct, gamma_bias, yerr = result_dic[key][3],
                ecolor='black', fmt='o', zorder=-500,markersize=1)
    # print(result_dic[key][3])
    plt.plot(np.linspace(0, len(result_dic)), np.linspace(0, len(result_dic))*0)
    plt.xlabel("ID",fontsize=27)
    plt.ylabel(which[1]+" bias (inferred - truth)",fontsize=27)
    plt.ylim(-0.4,0.4)
    plt.tick_params(labelsize=20)
    ct = ct+1
plt.show()
print(np.mean(gamma_bias_list), np.std(gamma_bias_list))
