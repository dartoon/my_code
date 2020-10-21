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
sys.path.insert(0,'../../../../py_tools/')
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

# # folder_type  = folder_list[0][:-3]
# # file_type = 'model_result_use_drz_Noisemap_subg3.pkl'
# file_type = 'model_result_calNoiseMap_modNoisemap_boostPossionx3_subg3.pkl'
# # file_type = 'model_result_calNoiseMap_modNoisemap_boostPossionx3_subg3.pkl'
file_type = 'result_modNoisemap_boostPossionx3_subg3.pkl'
folder_type = 'simulations_700_subg30/sim_lens_ID_subg30_'


# # # # file_type = 'model_result_subg3.pkl'
# # file_type = 'model_result_calNoiseMap_modNoisemap_boostPossionx8_noPSFerr_subg2_fixgamma.pkl'
# # # file_type = 'model_result_calNoiseMap_modNoisemap_useGrad_noPSFerr_subg3.pkl'
# file_type = 'result_calNoiseMap_modNoisemap_boostPossionx8_noPSFerr_subg3.pkl'
# folder_type = 'simulations_700_subg30/sim_lens_noqso_ID_subg30_'


folder_list = glob.glob(folder_type+'*')
folder_list.sort()
test_numer = 10 #len(folder_list)
folder_list = folder_list[:test_numer]

# for folder in folder_list:
# for folder in ['simulations_700_subg30/sim_lens_ID_subg30_718']:   #Gamma strongly biased.
# for folder in ['simulations_700_subg30/sim_lens_ID_subg30_724']:   #Chisq large
for folder in ['simulations_700_subg30/sim_lens_ID_subg30_758']: 
# for folder in ['simulations_700_subg30/sim_lens_ID_subg30_740']:  #larger gamma bias, but small H0 bias
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
    multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(folder+file_type,'rb'))
    # multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(folder+'model_result_use_drz_Noisemap_PSFnoisemapX0.1_subg3.pkl','rb'))

    lens_data = pyfits.getdata(folder+'Drz_QSO_image.fits')
    lens_mask = cr_mask(lens_data, 'normal_mask.reg')
    framesize = len(multi_band_list[0][0]['image_data'])  #81
    ct = int((len(lens_data) - framesize)/2)
    lens_data = lens_data[ct:-ct,ct:-ct]
    lens_mask = (1-lens_mask)[ct:-ct,ct:-ct]
    
    fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
    labels_new = [r"$\gamma$", r"$D_{\Delta t}$","H$_0$" ]    
    modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[lens_mask])
    f, axes = modelPlot.plot_main()
    plt.show()
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
    
    from lenstronomy.LensModel.lens_model import LensModel
    lens_model_list = ['PEMD','SHEAR']
    lens_model_class = LensModel(lens_model_list)    
    x_image, y_image = para_s[-1]['ra_image'], para_s[-1]['dec_image']
    source_pos = [para_s[2][0]['center_x'], para_s[2][0]['center_y']]
    potential = lens_model_class.potential(x_image, y_image, kwargs=kwargs_lens_list)
    fer_pot = lens_model_class.fermat_potential(x_image, y_image, kwargs_lens_list)
    geometry = potential + fer_pot

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
    for i in range(len(true_potential_diff)):
        print("potential mismatch:", true_potential_diff[i]- infe_potential_diff[i])
        print("geometry mismatch:", true_geometry_diff[i]- infe_geometry_diff[i])    
    # print("potential mismatch overall:", np.sqrt(np.sum(true_potential_diff- infe_potential_diff)**2) )
    # print("geometry mismatch overall :", np.sqrt(np.sum(true_geometry_diff- infe_geometry_diff)**2 ) )
    print("Fermat mismatch:", np.sum( (fer_pot[0] - fer_pot[1:])  - (fer_pot_infe[0] - fer_pot_infe[i]) ) ) 
    