#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 11:14:36 2020

@author: Dartoon
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
from lenstronomy.Plots.model_plot import ModelPlot
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.Util import constants as const
import scipy.optimize as op
from astropy.cosmology import FlatLambdaCDM
from lenstronomy.Cosmo.lens_cosmo import LensCosmo

file_type = 'model_result_noTD_subg3.pkl'
# file_type = 'model_result_calNoiseMap_modNoisemap_boostPossionx8_noPSFerr_subg3.pkl'
folder_type = 'simulations_700_subg30/sim_lens_noqso_ID_subg30_'

with_qso_savename = 'model_result_calNoiseMap_modNoisemap_boostPossionx3_subg3.pkl'

# folder_list = glob.glob(folder_type+'*')
# folder_list.sort()
# test_numer = 20
# folder_list = folder_list[1:test_numer]

folder_list = ['simulations_700_subg30/sim_lens_noqso_ID_subg30_702',
 # 'simulations_700_subg30/sim_lens_noqso_ID_subg30_703',
 # 'simulations_700_subg30/sim_lens_noqso_ID_subg30_704',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_705',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_706',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_708',
 # 'simulations_700_subg30/sim_lens_noqso_ID_subg30_709',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_710',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_712',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_713',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_714',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_715',
 # 'simulations_700_subg30/sim_lens_noqso_ID_subg30_716',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_717',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_718',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_720',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_721',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_722',
 'simulations_700_subg30/sim_lens_noqso_ID_subg30_724']

id_range = int(folder_list[0][-3:]), int(folder_list[-1][-3:])

def lnlike(theta, TD_obs, TD_err_l):
    tdd = theta
    TD = tdd/const.c * f_potential / const.day_s * const.arcsec**2 * const.Mpc
    TD = TD - TD[0]
    if tdd>0:
        return -0.5*( np.sum( (TD[1:]-TD_obs[1:])**2 / TD_err_l[1:]**2) ) 
    else:
        return -np.inf 
nll = lambda *args: -lnlike(*args)

def cal_Ddt(zl, zs, H0_ini=100, om=0.27):
    cosmo = FlatLambdaCDM(H0=H0_ini, Om0=0.27) 
    lensunits=LensCosmo(z_lens=zl, z_source=zs,cosmo= cosmo)
    D_l=lensunits.dd
    D_s=lensunits.ds
    D_ls=lensunits.dds
    Ddt_corr = (1+zl)*D_l*D_s/D_ls
    return Ddt_corr

def cal_h0(zl, zs, Ddt, om=0.27):
    Ddt_corr = cal_Ddt(zl, zs, H0_ini=100)
    ratio = Ddt_corr/Ddt
    return 100 * ratio
    
result_dic = {}
H0_list = []
for folder in folder_list:
    ID = folder[-3:]
    folder = folder+'/'
    read_file = folder+ file_type
    print(read_file)
    model_lists, para_s, lens_info= pickle.load(open(folder+'sim_kwargs.pkl','rb'))
    lens_model_list, lens_light_model_list, source_model_list, point_source_list = model_lists
    lens_model_list[0] = 'PEMD'
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
    multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(read_file,'rb'))
    fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
    kwargs_model['lens_model_list'][0] =  'PEMD'
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
    result_dic[folder+file_type] = [truth_dic, kwargs_result, chisq]
    #Take QSO postiion:
    qso_folder = folder_type[:-16] + 'ID_subg30_{0}/'.format(ID)
    _, _, kwargs_result_withQSO, _, _, _ = pickle.load(open(qso_folder+with_qso_savename,'rb'))
    lens_model_list = ['PEMD','SHEAR']
    lens_model_class = LensModel(lens_model_list)
    #Calculate the Fermat potential using QSO position and inferred len parameter
    # x_image, y_image = kwargs_result_withQSO['kwargs_ps'][0]['ra_image'], kwargs_result_withQSO['kwargs_ps'][0]['dec_image']  
    x_image, y_image = kwargs_ps['ra_image'], kwargs_ps['dec_image']
    f_potential = lens_model_class.fermat_potential(x_image, y_image, kwargs_lens=kwargs_result['kwargs_lens'])
    # TD = TD_distance/const.c * f_potential / const.day_s * const.arcsec**2 * const.Mpc
    result = op.minimize(nll, [TD_distance], args=(TD_obs, TD_err_l))   
    print(TD_distance, result["x"])
    H0 = cal_h0(z_l ,z_s, result["x"])
    H0_list.append(H0[0])

#%%
H0_true = 73.907
fig, ax = plt.subplots(figsize=(11,8))
for i in range(len(folder_list)):
    folder = folder_list[i]
    H0 = H0_list[i]
    ID = int(folder[-3:])
    key = folder_type + '{0}/'.format(ID) + file_type
    plt.scatter(ID, H0,
                c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
    # plt.text(ID, H0, repr(round(result_dic[key][-1], 1)),fontsize=15)
plt.plot(np.linspace(id_range[0]-1, id_range[1]+1), np.linspace(id_range[0]-1, id_range[1]+1)*0 + H0_true)
plt.xlabel("ID",fontsize=27)
plt.ylabel("$H_0$",fontsize=27)
plt.ylim(50,100)
plt.tick_params(labelsize=20)
plt.show()
del H0_list[0] 
print(np.mean(H0_list), np.std(H0_list))
# H0_list = np.array(H0_list)
# plt.hist(H0_list[:, 0])
# plt.show()    
        
#%%test parameter bias:
#para = 'theta_E'  #'gamma'
#para = 'gamma'
which = ['kwargs_lens', 'gamma']
#which = ['kwargs_source', 'center_y']
fig, ax = plt.subplots(figsize=(11,8))
gamma_bias_list = []
for folder in folder_list:
    ID = int(folder[-3:])
    key = folder_type + '{0}/'.format(ID) + file_type
    gamma_bias = result_dic[key][1][which[0]][0][which[1]] - result_dic[key][0][which[0]][0][which[1]]
    gamma_bias_list.append(gamma_bias)
    plt.scatter(ID, gamma_bias,
                c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
    # plt.errorbar(ct, gamma_bias, yerr = result_dic[key][3],
    #             ecolor='black', fmt='o', zorder=-500,markersize=1)
    # print(result_dic[key][3])
plt.plot(np.linspace(id_range[0]-1, id_range[1]+1), np.linspace(id_range[0]-1, id_range[1]+1)*0)
plt.xlabel("ID",fontsize=27)
plt.ylabel(which[1]+" bias (inferred - truth)",fontsize=27)
plt.ylim(-0.4,0.4)
plt.tick_params(labelsize=20)
plt.show()
print(np.mean(gamma_bias_list), np.std(gamma_bias_list))
