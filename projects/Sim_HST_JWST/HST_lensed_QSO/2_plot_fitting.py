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



folder_type = 'sim_lens_ID_'
idx = -1

folder_type = 'sim_lens_noqso_ID_'
idx = -2
     
#ID = 604
#ID = 605
#ID = 607
ID = 612

for ID in range(ID, ID+1):    
    folder = folder_type + '{0}/'.format(ID)
    files = glob.glob(folder+'model_result*.pkl')
    files.sort()
    read_file = files[idx]        
    
    folder = folder_type+'{0}/'.format(ID)
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

    #%%Print fitting result:
    labels_new = [r"$\gamma$", r"$D_{\Delta t}$","H$_0$" ]
    multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(read_file,'rb'))
    fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
#    labels_new = [r"$\gamma$", r"$D_{\Delta t}$","H$_0$" ]
    lens_mask = pyfits.getdata('noqso_mask.fits')
    modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat",likelihood_mask_list=[lens_mask])
    f, axes = modelPlot.plot_main()
    f.show()
#    f, axes = modelPlot.plot_separate()
#    f.show()
#    f, axes = modelPlot.plot_subtract_from_data_all()
#    f.show()
    plt.show()
##    for i in range(len(chain_list)):
##        chain_plot.plot_chain_list(chain_list, i)
##    plt.close()
    truths=[para_s[0][0]['gamma'],TD_distance, 73.907]	
    plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True, #range= [[0.8,1.5],[1,3],[0,1],[0, 1],[2000,5000],[20,100]], 
                         quantiles=[0.16, 0.5, 0.84], truths =truths,
                         title_kwargs={"fontsize": 15}, label_kwargs = {"fontsize": 25},
                         levels=1.0 - np.exp(-0.5 * np.array([1.,2.]) ** 2))
    plt.show()
    print('Truth:', truths)
