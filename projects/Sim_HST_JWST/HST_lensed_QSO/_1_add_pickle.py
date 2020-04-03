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

#folder_type = 'sim_lens_ID_'
folder_type = 'sim_lens_noqso_ID_'
result_dic = {}
id_range= [501, 522]
for ID in range(id_range[0], id_range[1]):  
    folder = folder_type + '{0}/'.format(ID)
    print(folder)
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
    if glob.glob(folder+'model_result.pkl') == []:
        result_dic[folder] = [None, None]
        continue
    multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting = pickle.load(open(folder+'2nd_model_result.pkl','rb'))
    fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
    
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.)
    from lenstronomy.Analysis.td_cosmography import TDCosmography
    td_cosmo = TDCosmography(z_l, z_s, kwargs_model, cosmo_fiducial=cosmo)
    from lenstronomy.Sampling.parameters import Param
    # make instance of parameter class with given model options, constraints and fixed parameters #
    param = Param(kwargs_model, fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo, 
                  kwargs_lens_init=kwargs_result['kwargs_lens'], **kwargs_constraints)
    sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[-1]
    mcmc_new_list = []
    labels_new = [r"$\gamma$", r"$D_{\Delta t}$","H$_0$" ]
    for i in range(len(samples_mcmc)):
        # transform the parameter position of the MCMC chain in a lenstronomy convention with keyword arguments #
        kwargs_result = param.args2kwargs(samples_mcmc[i])
        D_dt = kwargs_result['kwargs_special']['D_dt']
        fermat_pot = td_cosmo.fermat_potential(kwargs_result['kwargs_lens'], kwargs_result['kwargs_ps'])
    #    delta_fermat_12 = fermat_pot[0] - fermat_pot[2]
        gamma = kwargs_result['kwargs_lens'][0]['gamma']
    #    phi_ext, gamma_ext = kwargs_result['kwargs_lens'][1]['gamma1'], kwargs_result['kwargs_lens'][1]['gamma2']
        mcmc_new_list.append([gamma, D_dt, cal_h0(z_l ,z_s, D_dt)])        
    pickle.dump([multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list], open(folder+'2nd_model_result_newlist.pkl', 'wb'))
