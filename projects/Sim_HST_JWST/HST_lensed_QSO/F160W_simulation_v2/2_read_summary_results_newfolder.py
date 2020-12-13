#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 21:20:09 2017

@author: dxh
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
#file name:
filt='f160w'
import pickle
import sys
sys.path.insert(0,'../../../../py_tools/')
from flux_profile import cr_mask
from lenstronomy.Plots.model_plot import ModelPlot
from mask_objects import find_loc_max

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'


result_dic = {}


# file_type = 'result_PSFerr001_PSFinter_subg3.pkl'
# folder_type = 'simulations_700_subg30/sim_lens_ID_subg30_'


# file_type = 'result_PSFerr001_subg3.pkl'
# folder_type = 'simulations_700_subg30/sim_lens_noqso_ID_subg30_'


# folder_list = glob.glob(folder_type+'*')
# folder_list.sort()
# outlier = [736, 728, 748, 710, 715]
# folder_list = [folder_list[i] for i in range(len(folder_list)) if int(folder_list[i][-3:]) not in outlier]

# folder_type = 'AGN_result_folder/???_result_PSFinter_subg3.pkl'
folder_type = 'AGN_result_folder/idx*_ID*_centerPSF001_PSFinter.pkl'

folder_list = glob.glob(folder_type)
folder_list.sort()
# test_numer = len(folder_list)
# id_range = int(folder_list[0][-3:]), int(folder_list[-1][-3:])
id_range = [701, 770]

chisq_list = []
for folder in folder_list:
    folder = folder+'/'
    read_file = folder[:-1]
    ID = folder.split('ID')[1][:3]
    sim_folder = 'simulations_700_subg30/' + 'sim_lens_ID_subg30_' + ID +'/'
    print(read_file)
    model_lists, para_s, lens_info= pickle.load(open(sim_folder+'sim_kwargs.pkl','rb'))
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
    
    lens_data = pyfits.getdata(sim_folder+'Drz_QSO_image.fits')
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
    # for i in range(len(x_s)):
    #     lens_mask[np.sqrt((y_grid-y_s[i])**2 + (x_grid-x_s[i])**2) <4] = 0
        
    # plt.imshow(multi_band_list[0][0]['noise_map'],origin='lower')
    # plt.show()
    # plt.imshow(lens_mask,origin='lower')
    # plt.show()        
    modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[lens_mask])
    logL = modelPlot._imageModel.likelihood_data_given_model(source_marg=False, linear_prior=None, **kwargs_result)
    n_data = modelPlot._imageModel.num_data_evaluate
    chisq = -logL * 2 / n_data
    
    #!!!
    kwargs_result['kwargs_lens'][0]['gamma'] = np.median(chain_list[-1][1][:,0])
    
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
    chisq_list.append(chisq)

#%%
chisq_thre = 100 #np.percentile(chisq_list,80)
H0_true = 73.907
fig, ax = plt.subplots(figsize=(11,8))
H0_list = []
use_folder = []
for folder in folder_list:
    ID = folder.split('ID')[1][:3]
    key = folder
    # if abs(result_dic[key][-1]) < 3.0 and result_dic[key][2][1] < 90 and result_dic[key][2][1] >61:      #Use folders meets this requirments
    if abs(result_dic[key][-1]) < chisq_thre and result_dic[key][2][1] < 100 and result_dic[key][2][1] >60:      #Use folders meets this requirments
        ID = int(ID)
        H0 = result_dic[key][2]
        plt.scatter(ID, H0[1],
                    c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
        plt.errorbar(ID, H0[1], yerr = [[H0[2]-H0[1]], [H0[1]-H0[0]]],
                    ecolor='black', fmt='o', zorder=-500,markersize=1)  
        plt.text(ID, H0[1], repr(round(result_dic[key][-1], 3)),fontsize=15)
        H0_list.append([H0[1], (H0[2]-H0[1]+ H0[1]-H0[0])/2 ])
        use_folder.append(folder)
    else:
        print(ID)
#        ax.set_xticks(range(id_range[0]-1, id_range[1]+1,3)) 
plt.plot(np.linspace(id_range[0]-1, id_range[1]+1), np.linspace(id_range[0]-1, id_range[1]+1)*0 + H0_true, 'blue')
plt.xlabel("ID",fontsize=27)
plt.ylabel("$H_0$",fontsize=27)
plt.ylim(50,105)
plt.tick_params(labelsize=20)
plt.show()

H0_list = np.array(H0_list)
plt.hist(H0_list[:, 0])
plt.show()


submit_sum = np.array(H0_list)
print(np.mean(submit_sum[:,0]),'+-', np.std(submit_sum[:,0]))

# goodness_sum = round(1/float(len(submit_sum))*np.sum(((submit_sum[:,0]-H0_true)/submit_sum[:,1])**2),3)
# precision_sum = round(1/float(len(submit_sum))*np.sum(submit_sum[:,1]/H0_true)*100, 3)
# print("precision_sum", precision_sum)
# accuracy_sum = round(1/float(len(submit_sum))*np.sum((submit_sum[:,0]-H0_true)/H0_true)*100, 3)
# print("accuracy_sum", accuracy_sum)


#%%test parameter bias:
#para = 'theta_E'  #'gamma'
#para = 'gamma'
which = ['kwargs_lens', 'gamma']
#which = ['kwargs_source', 'center_y']
fig, ax = plt.subplots(figsize=(11,8))
gamma_bias_list = []
for folder in folder_list:
    ID = folder.split('ID')[1][:3]
    key = folder
    if abs(result_dic[key][-1]) < chisq_thre and result_dic[key][2][1] < 100 and result_dic[key][2][1] >60:      #Use folders meets this requirments
        gamma_bias = result_dic[key][1][which[0]][0][which[1]] - result_dic[key][0][which[0]][0][which[1]]
        gamma_bias_list.append(gamma_bias)
        plt.scatter(int(ID[1:]), gamma_bias,
                    c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
plt.plot(np.linspace(0, int(ID[1:])), np.linspace(0, int(ID[1:]))*0)
plt.xlabel("ID",fontsize=27)
plt.ylabel(which[1]+" bias (inferred - truth)",fontsize=27)
plt.ylim(-0.4,0.4)
plt.tick_params(labelsize=20)
plt.close()
print(which[1], np.mean(gamma_bias_list), np.std(gamma_bias_list))

#%%
if 'noqso' in folder_type:
    text = 'non-AGN'
else:
    text = 'AGN'
fig, ax = plt.subplots(figsize=(8,8))
for folder in folder_list:
    ID = folder.split('/')[1][:3]
    key = folder
    if abs(result_dic[key][-1]) < chisq_thre and result_dic[key][2][1] < 100 and result_dic[key][2][1] >60:      #Use folders meets this requirments
        H0 = result_dic[key][2]
        gamma_bias = result_dic[key][1][which[0]][0][which[1]] - result_dic[key][0][which[0]][0][which[1]]
        plt.scatter(gamma_bias, H0[1] - H0_true,
                    c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
        # plt.errorbar(gamma_bias, H0[1] - H0_true,
        #               yerr = [[H0[2]-H0[1]], [H0[1]-H0[0]]],
        #             ecolor='black', fmt='o', zorder=-500,markersize=1)      
        # plt.text(gamma_bias, H0[1] - H0_true, key[-3:],fontsize=15)
    #    ax.set_xticks(range(id_range[0]-1, id_range[1]+1,3)) 
plt.plot(np.linspace(-0.3,0.3), np.linspace(-0.3,0.3)*0, 'k', alpha = 0.5 )
plt.title('Lensed {0} case'.format(text), fontsize=27)    
plt.xlabel("$\gamma$ bias (inferred - truth)", fontsize=27)
plt.ylabel("$H_0$ bias (inferred - truth)", fontsize=27)
plt.ylim(-20,20)
plt.xlim(-0.7, 0.7)
# plt.tick_params(labelsize=30)
plt.tick_params(which='both', width=2, length = 7, labelsize=30)
# plt.savefig('bias_result_{0}.pdf'.format(text), bbox_inches = "tight")
plt.show()

# #%%test parameter bias on q. #!!! Not right because the current kwarg_result is not the best answer yet...
# which = ['kwargs_lens']
# #which = ['kwargs_source', 'center_y']
# fig, ax = plt.subplots(figsize=(11,8))
# ct = 0
# q_bias_list = []
# for folder in folder_list:
#     ID = folder[-3:]
#     key = folder_type + '{0}'.format(ID)
#     e1_true, e2_true = result_dic[key][0][which[0]][0]['e1'], result_dic[key][0][which[0]][0]['e2']
#     theta_true, q_true = param_util.ellipticity2phi_q(e1_true, e2_true)
    
#     e1_inf, e2_inf = result_dic[key][1][which[0]][0]['e1'], result_dic[key][1][which[0]][0]['e2']
#     theta_inf, q_inf = param_util.ellipticity2phi_q(e1_inf, e2_inf)
    
#     q_bias = q_inf - q_true
#     if abs(result_dic[key][-1]) < chisq_thre and result_dic[key][2][1] < 100 and result_dic[key][2][1] >60:      #Use folders meets this requirments
#         q_bias_list.append(q_bias)
#         plt.scatter(int(ID[1:]), q_bias,
#                     c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
# plt.plot(np.linspace(0, int(ID[1:])), np.linspace(0, int(ID[1:]))*0)
# plt.xlabel("ID",fontsize=27)
# plt.ylabel(" bias of q (inferred - truth)",fontsize=27)
# plt.tick_params(labelsize=20)
# plt.close()

# plt.hist(q_bias_list)
# plt.close()

# print(np.median(q_bias_list), np.std(q_bias_list))

# fig, ax = plt.subplots(figsize=(11,8))
# ct = 0
# for folder in folder_list:
#     ID = folder[-3:]
#     key = folder_type + '{0}'.format(ID)
#     H0 = result_dic[key][2]

#     e1_true, e2_true = result_dic[key][0][which[0]][0]['e1'], result_dic[key][0][which[0]][0]['e2']
#     theta_true, q_true = param_util.ellipticity2phi_q(e1_true, e2_true)
    
#     e1_inf, e2_inf = result_dic[key][1][which[0]][0]['e1'], result_dic[key][1][which[0]][0]['e2']
#     theta_inf, q_inf = param_util.ellipticity2phi_q(e1_inf, e2_inf)
    
#     q_bias = q_inf - q_true
#     if abs(result_dic[key][-1]) < chisq_thre and result_dic[key][2][1] < 100 and result_dic[key][2][1] >60:      #Use folders meets this requirments
#         plt.scatter(q_bias, H0[1] - H0_true,
#                     c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
#         # plt.errorbar(q_bias, H0[1] - H0_true,
#         #               yerr = [[H0[2]-H0[1]], [H0[1]-H0[0]]],
#         #             ecolor='black', fmt='o', zorder=-500,markersize=1)      
#         plt.text(q_bias, H0[1] - H0_true, key[-3:],fontsize=15)
#     #    ax.set_xticks(range(id_range[0]-1, id_range[1]+1,3)) 
# plt.xlabel("q bias (inferred - truth)", fontsize=27)
# plt.ylabel("$H_0$ bias (inferred - truth)", fontsize=27)
# plt.ylim(-20,30)
# plt.tick_params(labelsize=20)
# plt.show()

# #%%test parameter bias: #!!! Not right because the current kwarg_result is not the best answer yet...
# #para = 'theta_E'  #'gamma'
# #para = 'gamma'
# which = ['kwargs_lens', 'theta_E']
# #which = ['kwargs_source', 'center_y']
# fig, ax = plt.subplots(figsize=(11,8))
# param_bias_list = []
# for folder in folder_list:
#     ID = folder[-3:]
#     key = folder_type + '{0}'.format(ID)
#     if abs(result_dic[key][-1]) < chisq_thre and result_dic[key][2][1] < 100 and result_dic[key][2][1] >60:      #Use folders meets this requirments
#         param_bias = result_dic[key][1][which[0]][0][which[1]] - result_dic[key][0][which[0]][0][which[1]]
#         param_bias_list.append(param_bias)
#         plt.scatter(int(ID[1:]), param_bias,
#                     c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
# plt.plot(np.linspace(0, int(ID[1:])), np.linspace(0, int(ID[1:]))*0)
# plt.xlabel("ID",fontsize=27)
# plt.ylabel(which[1]+" bias (inferred - truth)",fontsize=27)
# plt.ylim(-0.4,0.4)
# plt.tick_params(labelsize=20)
# plt.close()
# print(which[1], np.mean(param_bias_list), np.std(param_bias_list))

# fig, ax = plt.subplots(figsize=(11,8))
# for folder in folder_list:
#     ID = folder[-3:]
#     key = folder_type + '{0}'.format(ID)
#     if abs(result_dic[key][-1]) < chisq_thre and result_dic[key][2][1] < 100 and result_dic[key][2][1] >60:      #Use folders meets this requirments
#         H0 = result_dic[key][2]
#         param_bias = result_dic[key][1][which[0]][0][which[1]] - result_dic[key][0][which[0]][0][which[1]]
#         plt.scatter(param_bias, H0[1] - H0_true,
#                     c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
#         # plt.errorbar(param_bias, H0[1] - H0_true,
#         #               yerr = [[H0[2]-H0[1]], [H0[1]-H0[0]]],
#         #             ecolor='black', fmt='o', zorder=-500,markersize=1)      
#         plt.text(param_bias, H0[1] - H0_true, key[-3:],fontsize=15)
#     #    ax.set_xticks(range(id_range[0]-1, id_range[1]+1,3)) 
# plt.xlabel(which[1]+" bias (inferred - truth)", fontsize=27)
# plt.ylabel("$H_0$ bias (inferred - truth)", fontsize=27)
# plt.ylim(-20,30)
# plt.tick_params(labelsize=20)
# plt.show()
