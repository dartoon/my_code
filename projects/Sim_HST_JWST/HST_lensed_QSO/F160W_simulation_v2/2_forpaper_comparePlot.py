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
import pickle
from lenstronomy.Plots.model_plot import ModelPlot
import sys
sys.path.insert(0,'../../../../py_tools/')
from flux_profile import cr_mask

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

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
        model_lists, para_s, lens_info= pickle.load(open(folder+'sim_kwargs.pkl','rb'))
        lens_model_list, lens_light_model_list, source_model_list, point_source_list = model_lists
        
        # multi_band_list, kwargs_model, kwargs_result_best, chain_list, fix_setting, mcmc_new_list = pickle.load(open(folder+file_type,'rb'))
        # fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
        # modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result_best,
        #                       arrow_size=0.02, cmap_string="gist_heat")
        if i == 0:
            run_file  = 'AGN_result_folder_run3times/idx11_ID714_PSFerr025_notPSFinter_morePSO_0.pkl'
        if i == 1:
            run_file = folder+file_type
        qso_folder = 'sim_lens_ID_{0}/'.format(ID)
        
        multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(run_file,'rb'))
        lens_data = pyfits.getdata(folder+'Drz_QSO_image.fits')
        len_std = pyfits.getdata(folder+'noise_map.fits')
        lens_mask = cr_mask(lens_data, 'normal_mask.reg')
        framesize = len(multi_band_list[0][0]['image_data'])  #81
        ct = int((len(lens_data) - framesize)/2)
        lens_data = lens_data[ct:-ct,ct:-ct]
        len_std = len_std[ct:-ct,ct:-ct]
        lens_mask = (1-lens_mask)[ct:-ct,ct:-ct]
        modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result,
                              arrow_size=0.02, cmap_string="gist_heat")
        model, error_map_, cov_param, param = modelPlot._imageModel.image_linear_solve(inv_bool=True, **kwargs_result)
        model = model[0]#[ct:-ct,ct:-ct]
        # error_map  = np.sqrt(len_std**2+np.abs(error_map_[0]))
        # ct = 2
        # lens_data = lens_data[ct:-ct,ct:-ct]
        # len_std = len_std[ct:-ct,ct:-ct]
        # lens_mask = lens_mask[ct:-ct,ct:-ct]
        # model = model[ct:-ct,ct:-ct]
        
        model_gene = modelPlot._band_plot_list[0]
        if i == 0:
            point_source_add = True
        else:
            point_source_add = False
        model_noarc = model_gene._bandmodel.image(model_gene._kwargs_lens_partial, model_gene._kwargs_source_partial,
                                                  model_gene._kwargs_lens_light_partial, model_gene._kwargs_ps_partial,
                                                   source_add=False, lens_light_add=True, point_source_add=point_source_add)
        # model_noarc = model_noarc[ct:-ct,ct:-ct]
        model_arc_res = lens_data - model_noarc

        if i == 0:
            lens_AGN_data, lens_AGN_std, lens_AGN_model, lens_AGN_arcres = lens_data, np.sqrt(len_std**2+np.abs(error_map_[0])), model, model_arc_res
        elif i ==1:
            lens_SN_data, lens_SN_std, lens_SN_model, lens_SN_arcres = lens_data, len_std, model, model_arc_res
    fig, axs = plt.subplots(2, 4, figsize=(15/3*3.5, 9/3*3.5))
    ((ax11, ax12, ax13, ax14), (ax21, ax22, ax23, ax24)) = axs
    
    f11 = ax11.imshow(lens_AGN_data, origin='lower', norm=LogNorm(), cmap = 'gist_heat', 
                    vmin=lens_AGN_data.min(), vmax=lens_AGN_data.max()*0.8)
    clb11 = fig.colorbar(f11, ax=ax11, shrink=0.55, pad=0.01, aspect=15) 
    clb11.ax.tick_params(labelsize=15) 
    ax11.set_title('Mock lensed AGN', fontsize = 20)
    ax11.set_xticks([])
    ax11.set_yticks([])
    
    f12 = ax12.imshow(lens_AGN_model, origin='lower', norm=LogNorm(), cmap = 'gist_heat', 
                    vmin=lens_AGN_data.min(), vmax=lens_AGN_data.max()*0.8)
    clb12 = fig.colorbar(f12, ax=ax12, shrink=0.55, pad=0.01, aspect=15) 
    clb12.ax.tick_params(labelsize=15) 
    ax12.set_title('Best-fit model', fontsize = 20)
    ax12.set_xticks([])
    ax12.set_yticks([])
    
    f13 = ax13.imshow(lens_AGN_arcres, origin='lower', norm=LogNorm(), cmap = 'gist_heat', 
                    vmin=lens_AGN_data.min(), vmax=lens_AGN_data.max()*0.8)
    clb13 = fig.colorbar(f13, ax=ax13, shrink=0.55, pad=0.01, aspect=15) 
    clb13.ax.tick_params(labelsize=15) 
    ax13.set_title('Deflector and AGNs subtracted', fontsize = 20)
    ax13.set_xticks([])
    ax13.set_yticks([])    
    
    f14 = ax14.imshow((lens_AGN_data-lens_AGN_model)/lens_AGN_std * lens_mask, origin='lower',
                    cmap = 'bwr', vmin=-6, vmax=6)
    clb14 = fig.colorbar(f14, ax=ax14, shrink=0.55, pad=0.01, aspect=15) 
    clb14.ax.tick_params(labelsize=15) 
    ax14.set_title('Normalized residuals', fontsize = 20)
    ax14.set_xticks([])
    ax14.set_yticks([])     
    
    
    f21 = ax21.imshow(lens_SN_data, origin='lower', norm=LogNorm(), cmap = 'gist_heat',
                    vmin=lens_AGN_data.min(), vmax=lens_AGN_data.max()*0.8)
    clb21 = fig.colorbar(f21, ax=ax21, shrink=0.55, pad=0.01, aspect=15) 
    clb21.ax.tick_params(labelsize=15) 
    ax21.set_title('Mock lensed transient', fontsize = 20)
    ax21.set_xticks([])
    ax21.set_yticks([])
    
    f22 = ax22.imshow(lens_SN_model, origin='lower', norm=LogNorm(), cmap = 'gist_heat',
                    vmin=lens_AGN_data.min(), vmax=lens_AGN_data.max()*0.8)
    clb22 = fig.colorbar(f22, ax=ax22, shrink=0.55, pad=0.01, aspect=15) 
    clb22.ax.tick_params(labelsize=15) 
    ax22.set_title('Best-fit model', fontsize = 20)
    ax22.set_xticks([])
    ax22.set_yticks([])
    
    f23 = ax23.imshow(lens_SN_arcres, origin='lower', norm=LogNorm(), cmap = 'gist_heat',
                    vmin=lens_AGN_data.min(), vmax=lens_AGN_data.max()*0.8)
    clb23 = fig.colorbar(f23, ax=ax23, shrink=0.55, pad=0.01, aspect=15) 
    clb23.ax.tick_params(labelsize=15) 
    ax23.set_title('Deflector subtracted', fontsize = 20)
    ax23.set_xticks([])
    ax23.set_yticks([])

    f24 = ax24.imshow((lens_SN_data-lens_SN_model)/lens_SN_std * lens_mask, origin='lower',
                    cmap = 'bwr', vmin=-6, vmax=6)
    clb24 = fig.colorbar(f24, ax=ax24, shrink=0.55, pad=0.01, aspect=15) 
    clb24.ax.tick_params(labelsize=15) 
    ax24.set_title('Normalized residuals', fontsize = 20)
    ax24.set_xticks([])
    ax24.set_yticks([])            
    
    # plt.tight_layout()
    plt.subplots_adjust(wspace=0.1, hspace=-0.3)
    plt.savefig('fitting_comparison.pdf')
    plt.show() 

# #%%
# test = modelPlot._band_plot_list[0]
# model = test._bandmodel.image(test._kwargs_lens_partial, test._kwargs_source_partial,
#                               test._kwargs_lens_light_partial,test._kwargs_ps_partial,
#                               lens_light_add=True, unconvolved=False, source_add = False)
# model = model[ct:-ct,ct:-ct]
# plt.imshow(lens_data-model, origin='low', norm=LogNorm(), cmap = 'gist_heat')
# plt.show()