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
    D_l=lensunits.dd
    D_s=lensunits.ds
    D_ls=lensunits.dds
    Ddt_corr = (1+zl)*D_l*D_s/D_ls
    return Ddt_corr

def cal_h0(zl, zs, Ddt, om=0.27):
    Ddt_corr = cal_Ddt(zl, zs, H0_ini=100)
    ratio = Ddt_corr/Ddt
    return 100 * ratio

test_numer = 50 #len(50)
kernel = 1
run_n = int(test_numer/kernel)

kernel_i = 0 # 0, 1 ,2, 3 .. max = kernel-1

# folder_list = glob.glob('simulations_700_subg30/sim_lens_ID_subg30_7??')
# folder_list.sort()
# savename = 'result_modNoisemap_boostPossionx3_subg3.pkl' #+ Simon's points; PSF not change, psf_error_map 0.1

folder_list = glob.glob('simulations_700_subg30/sim_lens_noqso_ID_subg30_7??')
folder_list.sort()
savename = 'result_calNoiseMap_modNoisemap_boostPossionx8_noPSFerr_subg3.pkl'


folder_list = folder_list[:test_numer]
#After talk with Simon:


# for folder in folder_list[kernel_i*run_n:kernel_i*run_n+run_n]:
for folder in [folder_list[-3]]:
    ID = folder[-3:]
    folder = folder + '/'
    print(folder)
    qso_folder = 'simulations_700_subg30/sim_lens_ID_subg30_{0}/'.format(ID)
    model_lists, para_s, lens_info= pickle.load(open(folder+'sim_kwargs.pkl','rb'))
    lens_model_list, lens_light_model_list, source_model_list, point_source_list = model_lists
    lens_model_list[0] = 'PEMD'
    z_l, z_s, TD_distance, TD_true, TD_obs, TD_err_l = lens_info
    kwargs_lens_list, kwargs_lens_light_list, kwargs_source_list, kwargs_ps = para_s
    solver_type = 'PROFILE_SHEAR'
    kwargs_constraints = {'joint_source_with_point_source': [[0, 0]],
                          'num_point_source_list': [len(kwargs_ps['ra_image'])],
                          'solver_type': solver_type,  # 'PROFILE', 'PROFILE_SHEAR', 'ELLIPSE', 'CENTER'
                          'Ddt_sampling': True,
                                  }
    lens_data = pyfits.getdata(folder+'Drz_QSO_image.fits')
    lens_mask = cr_mask(lens_data, 'normal_mask.reg')
    framesize = 81
    ct = int((len(lens_data) - framesize)/2)
    lens_data = lens_data[ct:-ct,ct:-ct]
    lens_mask = (1-lens_mask)[ct:-ct,ct:-ct]

    print(ID)
    multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(folder+savename,'rb'))
    fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
    labels_new = [r"$\gamma$", r"$D_{\Delta t}$","H$_0$" ]    
    modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat",
                          likelihood_mask_list=[lens_mask])
    f, axes = modelPlot.plot_main()
    f.show()
    # f, axes = modelPlot.plot_separate()
    # f.show()
    # f, axes = modelPlot.plot_subtract_from_data_all()
    # f.show()
    
    # for i in range(len(chain_list)):
    #     chain_plot.plot_chain_list(chain_list, i)
    # plt.show()
    
    # truths=[para_s[0][0]['gamma'],TD_distance, 73.907]	
    # plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True, #range= [[0.8,1.5],[1,3],[0,1],[0, 1],[2000,5000],[20,100]], 
    #                      quantiles=[0.16, 0.5, 0.84], truths =truths,
    #                      title_kwargs={"fontsize": 15}, label_kwargs = {"fontsize": 25},
    #                      levels=1.0 - np.exp(-0.5 * np.array([1.,2.]) ** 2))
    plt.show()
