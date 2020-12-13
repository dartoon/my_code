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

folder_type = 'AGN_result_folder/idx*_ID*_centerPSF001_PSFinter.pkl'

folder_list = glob.glob(folder_type)
folder_list.sort()


# index = int(sys.argv[1])
# print("Which index:", index)
savename = 'centerPSF001_PSFinter.pkl'

#%%plot_spec_filter
save_pkl_folder = 'AGN_result_folder/'
for folder in folder_list:
    ID = folder.split('ID')[1][:3]
    folder = folder + '/'
    print(folder)
    
    sim_folder = 'simulations_700_subg30/' + 'sim_lens_ID_subg30_' + ID +'/'
    model_lists, para_s, lens_info= pickle.load(open(sim_folder+'sim_kwargs.pkl','rb'))
    lens_model_list, lens_light_model_list, source_model_list, point_source_list = model_lists
    lens_model_list[0] = 'PEMD'
    z_l, z_s, TD_distance, TD_true, TD_obs, TD_err_l = lens_info
    kwargs_lens_list, kwargs_lens_light_list, kwargs_source_list, kwargs_ps = para_s
    solver_type = 'PROFILE_SHEAR'
    
    kwargs_constraints = {#'joint_source_with_point_source': [[0, 0]],
                          'num_point_source_list': [len(kwargs_ps['ra_image'])],
                          'solver_type': solver_type,  # 'PROFILE', 'PROFILE_SHEAR', 'ELLIPSE', 'CENTER'
                          'Ddt_sampling': True
                                  }
    index = folder.split('idx')[1].split('_')[0]
    save_file = save_pkl_folder+'idx{0}_ID'.format(index)+ID+'_'+savename
    multi_band_list, kwargs_model, kwargs_result_best, chain_list, fix_setting, mcmc_new_list = pickle.load(open(save_file,'rb'))
    # fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
    labels_new = [r"$\gamma$", r"$D_{\Delta t}$","H$_0$" ]
    modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result_best, arrow_size=0.02, cmap_string="gist_heat")
    f, axes = modelPlot.plot_main()
    plt.show()
    # f, axes = modelPlot.plot_separate()_
    # f.show()
    # f, axes = modelPlot.plot_subtract_from_data_all()
    # f.show()
    # multi_band_list = fitting_seq.multi_band_list
    # kwargs_psf_updated = multi_band_list[0][1]
    # f, axes = chain_plot.psf_iteration_compare(kwargs_psf_updated)
    # f.show()

#    for i in range(len(chain_list)):
#        chain_plot.plot_chain_list(chain_list, i)
#    plt.close()
    truths=[para_s[0][0]['gamma'],TD_distance, 73.907]	
    plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True, #range= [[0.8,1.5],[1,3],[0,1],[0, 1],[2000,5000],[20,100]], 
                          quantiles=[0.16, 0.5, 0.84], truths =truths,
                          title_kwargs={"fontsize": 15}, label_kwargs = {"fontsize": 25},
                          levels=1.0 - np.exp(-0.5 * np.array([1.,2.]) ** 2))
    plt.close()
