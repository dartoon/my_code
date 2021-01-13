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
#file name:kwargs_ps
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
folder_type = 'AGN_result_folder_run3times/idx*_ID*_PSFerr025_notPSFinter_morePSO_?.pkl'
# folder_type = 'AGN_result_folder_run3times/idx*_ID*noPSFerr_noPSFinter_morePSO*.pkl'


folder_list = glob.glob(folder_type)
folder_list.sort()
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
    modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, likelihood_mask_list=[lens_mask])
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

#%%
chisq_thre = 100 #np.percentile(chisq_list,80)
H0_true = 73.907
fig, ax = plt.subplots(figsize=(11,8))
H0_list = []
use_folder = []
for folder in folder_list:
    ID = folder.split('ID')[1][:3]
    key = folder
    runs = glob.glob(key.split('/')[0] + '/' + '*ID{0}_*'.format(ID) + folder_type[-33:])
    chisqs = [result_dic[runs[i]][-1] for i in range(len(runs))]
    best_chisq = np.min(chisqs)
    # if abs(result_dic[key][-1]) < chisq_thre and result_dic[key][2][1] < 90 and result_dic[key][2][1] >60:      #Use folders meets this requirments    
    if abs(result_dic[key][-1])  == best_chisq:      #Use folders meets this requirments    
        ID = int(ID)
        H0 = result_dic[key][2]
        plt.scatter(ID, H0[1],
                    c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
        plt.errorbar(ID, H0[1], yerr = [[H0[2]-H0[1]], [H0[1]-H0[0]]],
                    ecolor='black', fmt='o', zorder=-500,markersize=1)  
        # plt.text(ID, H0[1], repr(round(result_dic[key][-1], 3)),fontsize=15)
        H0_list.append([H0[1], (H0[2]-H0[1]+ H0[1]-H0[0])/2 ])
        use_folder.append(folder)
        
        
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

# #%%Test the accuracy and precision of the inferred point source
# acc, prec = [], []
# for folder_file in use_folder:
#     ID = folder_file.split('ID')[1][:3]
#     sim_folder = 'simulations_700_subg30/' + 'sim_lens_ID_subg30_' + ID +'/'
#     model_lists, para_s, lens_info= pickle.load(open(sim_folder+'sim_kwargs.pkl','rb'))
#     lens_model_list, lens_light_model_list, source_model_list, point_source_list = model_lists
#     z_l, z_s, TD_distance, TD_true, TD_obs, TD_err_l = lens_info
#     kwargs_lens_list, kwargs_lens_light_list, kwargs_source_list, kwargs_ps = para_s    
#     multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(folder_file,'rb'))
#     kwargs_ps_dec = [np.median(chain_list[-1][1][:,-5 + i]) for i in range(4)]        
#     kwargs_ps_dec = np.array(kwargs_ps_dec)
#     std_dec = np.std(kwargs_ps['dec_image'] - kwargs_ps_dec)  # Use the scatter of the general offset as the bias.  
#     kwargs_ps_ra = [np.median(chain_list[-1][1][:,-9 + i]) for i in range(4)]        
#     kwargs_ps_ra = np.array(kwargs_ps_ra)
#     std_ra= np.mean(kwargs_ps['dec_image'] - kwargs_ps_dec)  # Use the scatter of the general offset as the bias.
#     std_astrom = (std_dec ** 2 + std_ra ** 2 )** 0.5
#     acc.append(std_astrom)
#     mec_pre = [(np.std(chain_list[-1][1][:,-5 + i]) **2 + np.std(chain_list[-1][1][:,-9 +i]) **2) ** 0.5 for i in range(4)]
#     prec.append(mec_pre)

# print("accuracy and precision", np.mean(acc), np.mean(prec)) #accuracy is actually the std of the PS offset.

#%%test parameter bias:
#para = 'theta_E'  #'gamma'
#para = 'gamma'
which = ['kwargs_lens', 'gamma']
# #which = ['kwargs_source', 'center_y']
# fig, ax = plt.subplots(figsize=(11,8))
# gamma_bias_list = []
# for folder in folder_list:
#     ID = folder.split('ID')[1][:3]
#     key = folder
#     # if abs(result_dic[key][-1]) < chisq_thre and result_dic[key][2][1] < 100 and result_dic[key][2][1] >60:      #Use folders meets this requirments
#     runs = glob.glob(key.split('/')[0] + '/' + '*ID{0}*'.format(ID) + folder_type[-33:])
#     chisqs = [result_dic[runs[i]][-1] for i in range(len(runs))]
#     best_chisq = np.min(chisqs)
#     # if abs(result_dic[key][-1]) < chisq_thre and result_dic[key][2][1] < 90 and result_dic[key][2][1] >60:      #Use folders meets this requirments    
#     if result_dic[key][-1]  == best_chisq:      #Use folders meets this requirments  
#         gamma_bias = result_dic[key][1][which[0]][0][which[1]] - result_dic[key][0][which[0]][0][which[1]]
#         gamma_bias_list.append(gamma_bias)
#         plt.scatter(int(ID), gamma_bias,
#                     c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
# plt.plot(np.linspace(id_range[0]-1, id_range[1]+1), np.linspace(id_range[0]-1, id_range[1]+1)*0, 'blue')
# plt.xlabel("ID",fontsize=27)
# plt.ylabel(which[1]+" bias (inferred - truth)",fontsize=27)
# plt.ylim(-0.4,0.4)
# plt.tick_params(labelsize=20)
# plt.show()
# print(which[1], np.mean(gamma_bias_list), np.std(gamma_bias_list))

if 'noqso' in folder_type:
    text = 'non-AGN'
else:
    text = 'AGN'
fig, ax = plt.subplots(figsize=(8,8))
for folder in folder_list:
    ID = folder.split('ID')[1][:3]
    key = folder
    runs = glob.glob(key.split('/')[0] + '/' + '*ID{0}*'.format(ID) + folder_type[-33:])
    chisqs = [result_dic[runs[i]][-1] for i in range(len(runs))]
    best_chisq = np.min(chisqs)
    if abs(result_dic[key][-1])  == best_chisq:      #Use folders meets this requirments          H0 = result_dic[key][2]
        H0 = result_dic[key][2]        
        gamma_bias = result_dic[key][1][which[0]][0][which[1]] - result_dic[key][0][which[0]][0][which[1]]
        plt.scatter(gamma_bias, H0[1] - H0_true,
                    c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
plt.plot(np.linspace(-0.3,0.3), np.linspace(-0.3,0.3)*0, 'k', alpha = 0.5 )
plt.title('Lensed {0} case'.format(text), fontsize=27)    
plt.xlabel("$\gamma$ bias (inferred - truth)", fontsize=27)
plt.ylabel("$H_0$ bias (inferred - truth)", fontsize=27)
plt.ylim(-20,20)
plt.xlim(-0.25, 0.25)
# plt.tick_params(labelsize=30)
plt.tick_params(which='both', width=2, length = 7, labelsize=30)
# plt.savefig('bias_result_{0}.pdf'.format(text), bbox_inches = "tight")
plt.show()

import lenstronomy.Util.param_util as param_util
# #%%test parameter bias on q. #!!! Not right because the current kwarg_result is not the best answer yet...
which = ['kwargs_lens']
# #which = ['kwargs_source', 'center_y']
# fig, ax = plt.subplots(figsize=(11,8))
# ct = 0
# q_bias_list = []
# for folder in folder_list:
#     ID = folder.split('ID')[1][:3]
#     key = folder
#     e1_true, e2_true = result_dic[key][0][which[0]][0]['e1'], result_dic[key][0][which[0]][0]['e2']
#     theta_true, q_true = param_util.ellipticity2phi_q(e1_true, e2_true)
    
#     e1_inf, e2_inf = result_dic[key][1][which[0]][0]['e1'], result_dic[key][1][which[0]][0]['e2']
#     theta_inf, q_inf = param_util.ellipticity2phi_q(e1_inf, e2_inf)
    
#     q_bias = q_inf - q_true
#     if abs(result_dic[key][-1]) < chisq_thre and result_dic[key][2][1] < 100 and result_dic[key][2][1] >60:      #Use folders meets this requirments
#         q_bias_list.append(q_bias)
#         plt.scatter(int(ID), q_bias,
#                     c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
# plt.plot(np.linspace(id_range[0]-1, id_range[1]+1), np.linspace(id_range[0]-1, id_range[1]+1)*0, 'blue')
# plt.xlabel("ID",fontsize=27)
# plt.ylabel(" bias of q (inferred - truth)",fontsize=27)
# plt.tick_params(labelsize=20)
# plt.show()
# print(np.mean(q_bias_list), np.std(q_bias_list))
q_true_list = []
if 'noqso' in folder_type:
    text = 'non-AGN'
else:
    text = 'AGN'
fig, ax = plt.subplots(figsize=(8,8))
for folder in folder_list:
    ID = folder.split('ID')[1][:3]
    key = folder
    runs = glob.glob(key.split('/')[0] + '/' + '*ID{0}*'.format(ID) + folder_type[-33:])
    chisqs = [result_dic[runs[i]][-1] for i in range(len(runs))]
    best_chisq = np.min(chisqs)
    if abs(result_dic[key][-1])  == best_chisq:      #Use folders meets this requirments          H0 = result_dic[key][2]
        H0 = result_dic[key][2]        
        e1_true, e2_true = result_dic[key][0][which[0]][0]['e1'], result_dic[key][0][which[0]][0]['e2']
        theta_true, q_true = param_util.ellipticity2phi_q(e1_true, e2_true)
        
        e1_inf, e2_inf = result_dic[key][1][which[0]][0]['e1'], result_dic[key][1][which[0]][0]['e2']
        theta_inf, q_inf = param_util.ellipticity2phi_q(e1_inf, e2_inf)
        q_true_list.append(q_true)
        q_bias = q_inf - q_true
        plt.scatter(q_bias, H0[1] - H0_true,
                    c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
plt.plot(np.linspace(-0.3,0.3), np.linspace(-0.3,0.3)*0, 'k', alpha = 0.5 )
plt.title('Lensed {0} case'.format(text), fontsize=27)    
plt.xlabel("$q$ bias (inferred - truth)", fontsize=27)
plt.ylabel("$H_0$ bias (inferred - truth)", fontsize=27)
plt.ylim(-20,20)
plt.xlim(-0.08, 0.08)
# plt.tick_params(labelsize=30)
plt.tick_params(which='both', width=2, length = 7, labelsize=30)
# plt.savefig('bias_result_{0}.pdf'.format(text), bbox_inches = "tight")
plt.show()

#%%
if 'noqso' in folder_type:
    text = 'non-AGN'
else:
    text = 'AGN'
fig, ax = plt.subplots(figsize=(8,8))
for folder in folder_list:
    ID = folder.split('ID')[1][:3]
    key = folder
    runs = glob.glob(key.split('/')[0] + '/' + '*ID{0}*'.format(ID) + folder_type[-33:])
    chisqs = [result_dic[runs[i]][-1] for i in range(len(runs))]
    best_chisq = np.min(chisqs)
    if abs(result_dic[key][-1])  == best_chisq:      #Use folders meets this requirments          H0 = result_dic[key][2]
        H0 = result_dic[key][2]        
        e1_true, e2_true = result_dic[key][0][which[0]][0]['e1'], result_dic[key][0][which[0]][0]['e2']
        theta_true, q_true = param_util.ellipticity2phi_q(e1_true, e2_true)
        
        e1_inf, e2_inf = result_dic[key][1][which[0]][0]['e1'], result_dic[key][1][which[0]][0]['e2']
        theta_inf, q_inf = param_util.ellipticity2phi_q(e1_inf, e2_inf)
        q_bias = q_inf - q_true
        plt.scatter(q_true, H0[1] - H0_true,
                    c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
plt.plot(np.linspace(-0.3,1), np.linspace(-0.3,0.3)*0, 'k', alpha = 0.5 )
plt.title('Lensed {0} case'.format(text), fontsize=27)    
plt.xlabel("$q$ value of mass", fontsize=27)
plt.ylabel("$H_0$ bias (inferred - truth)", fontsize=27)
plt.ylim(-20,20)
plt.xlim(0.6, 1)
# plt.tick_params(labelsize=30)
plt.tick_params(which='both', width=2, length = 7, labelsize=30)
# plt.savefig('bias_result_{0}.pdf'.format(text), bbox_inches = "tight")
plt.show()

#Anowar's 2020 q distribution: plt.hist([0.69, 0.65, 0.89, 0.69, 0.59, 0.68, 0.72, 0.71, 0.56, 0.82, 0.91, 0.66, 0.64, 0.70, 0.68, 0.77, 0.87, 0.83, 0.66, 0.70, 0.69, 0.59, 0.65])