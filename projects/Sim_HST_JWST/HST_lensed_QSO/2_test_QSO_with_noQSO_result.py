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
import pickle
import corner    
from astropy.cosmology import FlatLambdaCDM
from lenstronomy.Workflow.fitting_sequence import FittingSequence
from lenstronomy.Cosmo.lens_cosmo import LensCosmo
from lenstronomy.Plots import chain_plot
from lenstronomy.Plots.model_plot import ModelPlot
from lenstronomy.Analysis.td_cosmography import TDCosmography
from lenstronomy.Sampling.parameters import Param
from lenstronomy.Data.psf import PSF
import lenstronomy.Util.simulation_util as sim_util
from lenstronomy.Data.imaging_data import ImageData
import sys
sys.path.insert(0,'../../../py_tools/')
from flux_profile import cr_mask
from mask_objects import find_loc_max
#file name:
filt='f160w'

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

# folder_list = glob.glob('sim_lens_noqso_ID_7??')
# folder_list.sort()
# test_numer = 1
# kernel = 1
# run_n = int(test_numer/kernel)

# kernel_i = 0 # 0, 1 ,2, 3 .. max = kernel-1

# folder_list = folder_list[1:2]
folder_list = ['simulations_700_subg30/sim_lens_ID_subg30_724']
# folder_list = ['simulations_700_subg30/sim_lens_noqso_ID_subg30_724']


for folder in folder_list:
    ID = folder[-3:]
    folder = folder + '/'
    print(folder)
    qso_folder = 'sim_lens_ID_{0}/'.format(ID)
    model_lists, para_s, lens_info= pickle.load(open(folder+'sim_kwargs.pkl','rb'))
    lens_model_list, lens_light_model_list, source_model_list, point_source_list = model_lists
    # lens_model_list[0] = 'PEMD'
    z_l, z_s, TD_distance, TD_true, TD_obs, TD_err_l = lens_info
    kwargs_lens_list, kwargs_lens_light_list, kwargs_source_list, kwargs_ps = para_s
    solver_type = 'PROFILE_SHEAR'
    
#    if len(kwargs_ps['ra_image']) <4:  #Would delete all the double
#        print(folder)
#        import shutil
#        shutil.rmtree(folder)
##        continue
##        if abs(kwargs_ps['ra_image']).min() != abs(kwargs_ps['ra_image'][-1]) and abs(kwargs_ps['dec_image']).min() != abs(kwargs_ps['dec_image'][-1]):
##            raise ValueError("The double image is not taken the points position correctly")
##        kwargs_ps['ra_image'] = kwargs_ps['ra_image'][:2] 
##        kwargs_ps['dec_image'] = kwargs_ps['dec_image'][:2]
##        kwargs_ps['point_amp'] = kwargs_ps['point_amp'][:2]
##        TD_obs = TD_obs[:2]
##        TD_err_l = TD_err_l[:2]
##        solver_type = 'THETA_E_PHI'
    kwargs_constraints = {'joint_source_with_point_source': [[0, 0]],
                          'num_point_source_list': [len(kwargs_ps['ra_image'])],
                          'solver_type': solver_type,  # 'PROFILE', 'PROFILE_SHEAR', 'ELLIPSE', 'CENTER'
                          'Ddt_sampling': True,
                                  }
    multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(folder+'model_result_use_drz_Noisemap_subg3.pkl','rb'))
    multi_band_list_noqso, kwargs_model_noqso, kwargs_result_noqso, chain_list_noqso, fix_setting_noqso, mcmc_new_list_noqso = pickle.load(open('simulations_700_subg30/sim_lens_noqso_ID_subg30_724/'+'model_result_use_drz_Noisemap_subg3.pkl','rb'))
    kwargs_result_noqso['kwargs_ps'] = kwargs_result['kwargs_ps']
    kwargs_result = kwargs_result_noqso

    fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
    labels_new = [r"$\gamma$", r"$D_{\Delta t}$","H$_0$" ]    
    modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat")
    f, axes = modelPlot.plot_main()
    f.show()
    # f, axes = modelPlot.plot_separate()
    # f.show()
    # f, axes = modelPlot.plot_subtract_from_data_all()
    # f.show()
    
    # for i in range(len(chain_list)):
    #     chain_plot.plot_chain_list(chain_list, i)
    # plt.show()
    
    truths=[para_s[0][0]['gamma'],TD_distance, 73.907]	
    plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True, #range= [[0.8,1.5],[1,3],[0,1],[0, 1],[2000,5000],[20,100]], 
                          quantiles=[0.16, 0.5, 0.84], truths =truths,
                          title_kwargs={"fontsize": 15}, label_kwargs = {"fontsize": 25},
                          levels=1.0 - np.exp(-0.5 * np.array([1.,2.]) ** 2))
    plt.show()

    # sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[-1]
    # param = Param(kwargs_model, fixed_lens, fixed_source, fixed_lens_light,
    #               kwargs_lens_init=kwargs_result['kwargs_lens'], **kwargs_constraints)        
    # mcmc_new_list = []     
    # steps = len(samples_mcmc)
    # for i in range(steps-5000, steps):
    #     kwargs_result = param.args2kwargs(samples_mcmc[i])
    #     thetaE = kwargs_result['kwargs_lens'][0]['theta_E']
    #     gamma = kwargs_result['kwargs_lens'][0]['gamma']
    #     e1 = kwargs_result['kwargs_lens'][0]['e1']
    #     e2 = kwargs_result['kwargs_lens'][0]['e2']
    #     mcmc_new_list.append([thetaE, gamma, e1, e2])    
    # plt.show()
    # labels_new = ["Theta_E", r"$\gamma$", r"$e1$", r"$e2$"] 
    # truths=[para_s[0][0]['theta_E'], para_s[0][0]['gamma'],para_s[0][0]['e1'], para_s[0][0]['e2']]	    
    # plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True, #range= [[0.8,1.5],[1,3],[0,1],[0, 1],[2000,5000],[20,100]], 
    #                       quantiles=[0.16, 0.5, 0.84], truths =truths,
    #                       title_kwargs={"fontsize": 15}, label_kwargs = {"fontsize": 25},
    #                       levels=1.0 - np.exp(-0.5 * np.array([1.,2.]) ** 2))
    # plt.show()
    