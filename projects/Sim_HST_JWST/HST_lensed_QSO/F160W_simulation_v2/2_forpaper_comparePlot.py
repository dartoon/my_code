#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 11:34:34 2020

@author: Xuheng Ding
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
# ID = str(707)

file_type_l = ['result_PSFerr001_PSFinter_subg3.pkl', 'result_PSFerr001_subg3.pkl']
folder_type_l = ['simulations_700_subg30/sim_lens_ID_subg30_', 'simulations_700_subg30/sim_lens_noqso_ID_subg30_']


# ID_list = [702 ,703 ,705 ,706 ,707 ,708 ,709 ,710 ,711 ,712 ,713, 714 ,715 ,717 ,718 ,719 ,720 ,721 ,724 ,727 ,728 ,729 ,730 ,731 ,733 ,734 ,735 ,736 ,737 ,738 ,739 ,740 ,741 ,742 ,743 ,744 ,745 ,747 ,748 ,749 ,751 ,752 ,753 ,754 ,755 ,757 ,759 ,760]

ID_list = [714]
for ID in ID_list:
    ID = str(ID)
    print(ID)
    for i in range(len(file_type_l)):
        file_type, folder_type = file_type_l[i], folder_type_l[i]
        folder = folder_type+ ID + '/'
        qso_folder = 'sim_lens_ID_{0}/'.format(ID)
        model_lists, para_s, lens_info= pickle.load(open(folder+'sim_kwargs.pkl','rb'))
        lens_model_list, lens_light_model_list, source_model_list, point_source_list = model_lists
        multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(folder+file_type,'rb'))
        lens_data = pyfits.getdata(folder+'Drz_QSO_image.fits')
        len_std = pyfits.getdata(folder+'noise_map.fits')
        lens_mask = cr_mask(lens_data, 'normal_mask.reg')
        framesize = len(multi_band_list[0][0]['image_data'])  #81
        ct = int((len(lens_data) - framesize)/2)
        lens_data = lens_data[ct:-ct,ct:-ct]
        len_std = len_std[ct:-ct,ct:-ct]
        lens_mask = (1-lens_mask)[ct:-ct,ct:-ct]
        modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[lens_mask])
        model, error_map, cov_param, param = modelPlot._imageModel.image_linear_solve(inv_bool=True, **kwargs_result)
        model = model[0]#[ct:-ct,ct:-ct]
        ct = 2
        lens_data = lens_data[ct:-ct,ct:-ct]
        len_std = len_std[ct:-ct,ct:-ct]
        lens_mask = lens_mask[ct:-ct,ct:-ct]
        model = model[ct:-ct,ct:-ct]
        
        if i == 0:
            lens_AGN_data, lens_AGN_std, lens_AGN_model = lens_data, len_std, model
        elif i ==1:
            lens_SN_data, lens_SN_std, lens_SN_model = lens_data, len_std, model
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    ((ax1, ax2), (ax3, ax4)) = axs
    
    f1 = ax1.imshow(lens_AGN_data, origin='lower', norm=LogNorm(), cmap = 'gist_heat', 
                    vmin=lens_AGN_data.min(), vmax=lens_AGN_data.max()*0.8)
    clb1 = fig.colorbar(f1, ax=ax1, shrink=0.9, pad=0.01, aspect=15) 
    clb1.ax.tick_params(labelsize=15) 
    ax1.set_title('Lensed AGN case', fontsize = 20)
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    f2 = ax2.imshow((lens_AGN_data-lens_AGN_model)/lens_AGN_std * lens_mask, origin='lower',
                    cmap = 'bwr', vmin=-6, vmax=6)
    clb2 = fig.colorbar(f2, ax=ax2, shrink=0.9, pad=0.01, aspect=15) 
    clb2.ax.tick_params(labelsize=15) 
    ax2.set_title('Residuals map', fontsize = 20)
    ax2.set_xticks([])
    ax2.set_yticks([])
    
    f3 = ax3.imshow(lens_SN_data, origin='lower', norm=LogNorm(), cmap = 'gist_heat',
                    vmin=lens_AGN_data.min(), vmax=lens_AGN_data.max()*0.8)
    clb3 = fig.colorbar(f3, ax=ax3, shrink=0.9, pad=0.01, aspect=15) 
    clb3.ax.tick_params(labelsize=15) 
    ax3.set_title('Lensed SNe case', fontsize = 20)
    ax3.set_xticks([])
    ax3.set_yticks([])
    
    f4 = ax4.imshow((lens_SN_data-lens_SN_model)/lens_SN_std * lens_mask, origin='lower',
                    cmap = 'bwr', vmin=-6, vmax=6)
    clb4 = fig.colorbar(f4, ax=ax4, shrink=0.9, pad=0.01, aspect=15) 
    clb4.ax.tick_params(labelsize=15) 
    ax4.set_title('Residuals map', fontsize = 20)
    ax4.set_xticks([])
    ax4.set_yticks([])
    
    plt.tight_layout()
    plt.savefig('fitting_comparison.pdf')
    plt.show() 
