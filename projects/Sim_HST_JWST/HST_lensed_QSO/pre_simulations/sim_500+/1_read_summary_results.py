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
from mask_objects import find_loc_max

from lenstronomy.Cosmo.lens_cosmo import LensCosmo
def cal_Ddt(zl, zs, H0_ini=100, om=0.27):
    cosmo = FlatLambdaCDM(H0=H0_ini, Om0=0.27) 
    lensunits=LensCosmo(z_lens=zl, z_source=zs,cosmo= cosmo)
    D_l=lensunits.D_d
    D_s=lensunits.D_s
    D_ls=lensunits.D_ds 
    Ddt_corr = (1+zl)*D_l*D_s/D_ls
    return Ddt_corr

def cal_h0(zl, zs, Ddt, om=0.27):
    Ddt_corr = cal_Ddt(zl, zs, H0_ini=100)
    ratio = Ddt_corr/Ddt
    return 100 * ratio

result_dic = {}
#folder_type = 'sim_lens_ID_'
#file_type = '2nd_model_result_newlist.pkl'

#folder_type = 'sim_lens_noqso_ID_'
#file_type = '2nd_model_result_improve.pkl'
id_range= [501, 522]
for ID in range(id_range[0], id_range[1]):  
    folder = folder_type + '{0}/'.format(ID)
#    print(folder)
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
    if glob.glob(folder+'2nd_model_result_newlist.pkl') == []:
        result_dic[folder] = [None, None]
        continue
    multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(folder+file_type,'rb'))
    fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
    mcmc_new_list = np.array(mcmc_new_list)
    H0_list = mcmc_new_list[:,2]
    truth_dic = {}
    truth_dic['kwargs_lens'] =kwargs_lens_list
    truth_dic['kwargs_source'] =kwargs_source_list
    truth_dic['kwargs_lens_light'] =kwargs_lens_light_list
    truth_dic['kwargs_ps'] = kwargs_ps
    truth_dic['D_dt'] = TD_distance
    result_dic[folder[:-1]] = [truth_dic, kwargs_result, [np.percentile(H0_list,16), np.percentile(H0_list,50), np.percentile(H0_list,84)]]

#%%
H0_true = 70.656 
fig, ax = plt.subplots(figsize=(11,8))
for ID in range(id_range[0], id_range[1]):
    key = folder_type + '{0}'.format(ID)
    H0 = result_dic[key][-1]
    plt.scatter(ID, H0[1],
                c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
    ax.set_xticks(range(id_range[0]-1, id_range[1]+1,3)) 
    plt.plot(np.linspace(id_range[0]-1, id_range[1]+1), np.linspace(id_range[0]-1, id_range[1]+1)*0 + H0_true)
    plt.xlabel("ID",fontsize=27)
    plt.ylabel("$H_0$",fontsize=27)
    plt.ylim(25,90)
    plt.tick_params(labelsize=20)
plt.show()

#%%
fig, ax = plt.subplots(figsize=(11,8))
for ID in range(id_range[0], id_range[1]):
    key = folder_type + '{0}'.format(ID)
    gamma_bias = result_dic[key][1]['kwargs_lens'][0]['gamma'] - result_dic[key][0]['kwargs_lens'][0]['gamma']
    plt.scatter(ID, gamma_bias,
                c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
    ax.set_xticks(range(id_range[0]-1, id_range[1]+1,3)) 
    plt.plot(np.linspace(id_range[0]-1, id_range[1]+1), np.linspace(id_range[0]-1, id_range[1]+1)*0)
    plt.xlabel("ID",fontsize=27)
    plt.ylabel("$\gamma$ bias (inferred - truth)",fontsize=27)
    plt.ylim(-0.23,0.23)
    plt.tick_params(labelsize=20)
plt.show()
