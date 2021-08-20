#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 21:20:09 2017

@author: dxh
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
#file name:
filt='f160w'
import pickle
import sys
sys.path.insert(0,'../../../../py_tools/')
from flux_profile import cr_mask
from lenstronomy.Plots.model_plot import ModelPlot
from mask_objects import find_loc_max

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'


result_dic = {}


# file_type = 'result_PSFerr001_PSFinter_subg3.pkl'
# folder_type = 'simulations_700_subg30/sim_lens_ID_subg30_'


file_type = 'result_PSFerr001_subg3.pkl'
folder_type = 'simulations_700_subg30/sim_lens_noqso_ID_subg30_'


folder_list = glob.glob(folder_type+'*')
folder_list.sort()
# outlier = [736, 728, 748, 710, 715]
# folder_list = [folder_list[i] for i in range(len(folder_list)) if int(folder_list[i][-3:]) not in outlier]


test_numer = len(folder_list)
folder_list = folder_list[:test_numer]

id_range = int(folder_list[0][-3:]), int(folder_list[-1][-3:])

chisq_list = []
for folder in folder_list:
    folder = folder+'/'
    read_file = folder+ file_type
    print(read_file)
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
#    if glob.glob(folder+'model_result.pkl') == []:
#        result_dic[folder] = [None, None]
#        continue
    multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(read_file,'rb'))
    fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
    mcmc_new_list = np.array(mcmc_new_list)
    H0_list = mcmc_new_list[:,-1]
    
    lens_data = pyfits.getdata(folder+'Drz_QSO_image.fits')
    lens_mask = cr_mask(lens_data, 'normal_mask.reg')
    framesize = len(multi_band_list[0][0]['image_data'])  #81
    ct = int((len(lens_data) - framesize)/2)
    lens_data = lens_data[ct:-ct,ct:-ct]
    lens_mask = (1-lens_mask)[ct:-ct,ct:-ct]
    
    x, y =find_loc_max(lens_data)
    x_s, y_s = [], []
    center = (framesize-1)/2
    deltaPix = 0.08
    for i in range(len(x)):
        x0, y0 =  (float(x[i]) - center) * deltaPix , (float(y[i]) - center) * deltaPix
#            print(x0, y0)
        ds = (x0- kwargs_ps['ra_image'] )**2 + (y0-kwargs_ps['dec_image'])**2
        if ds.min()<0.01:
            x_s.append(x[i])
            y_s.append(y[i])
    y_grid,x_grid = np.indices((framesize,framesize))   #with higher resolution 60*6
    # for i in range(len(x_s)):
    #     lens_mask[np.sqrt((y_grid-y_s[i])**2 + (x_grid-x_s[i])**2) <4] = 0
        
    # plt.imshow(multi_band_list[0][0]['noise_map'],origin='lower')
    # plt.show()
    # plt.imshow(lens_mask,origin='lower')
    # plt.show()        
    modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[lens_mask])
    logL = modelPlot._imageModel.likelihood_data_given_model(source_marg=False, linear_prior=None, **kwargs_result)
    n_data = modelPlot._imageModel.num_data_evaluate
    chisq = -logL * 2 / n_data
    
    #!!!
    kwargs_result['kwargs_lens'][0]['gamma'] = np.median(chain_list[-1][1][:,0])
    
    if chisq< 0.0:
        print(folder[-4:-1], round(np.median(H0_list), 3))
        f, axes = modelPlot.plot_main()
        f.show()
        plt.show()
    truth_dic = {}
    truth_dic['kwargs_lens'] =kwargs_lens_list
    truth_dic['kwargs_source'] =kwargs_source_list
    truth_dic['kwargs_lens_light'] =kwargs_lens_light_list
    truth_dic['kwargs_ps'] = kwargs_ps
    truth_dic['D_dt'] = TD_distance
    result_dic[folder[:-1]] = [truth_dic, kwargs_result, 
                               [np.percentile(H0_list,16), np.percentile(H0_list,50), np.percentile(H0_list,84)], 
                               (np.percentile(chain_list[-1][1][:,0],84) - np.percentile(chain_list[-1][1][:,0],16))/2, chisq]
    chisq_list.append(chisq)

#%%
import lenstronomy.Util.simulation_util as sim_util
import lenstronomy.Util.image_util as image_util
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.Data.psf import PSF
from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.ImSim.image_model import ImageModel
# data specifics
numPix = 201  # cutout pixel size
kwargs_data = sim_util.data_configure_simple(numPix, deltaPix, 100, .05)
data_class = ImageData(**kwargs_data)
kwargs_psf = {'psf_type': 'GAUSSIAN', 'fwhm': 0.16, 'pixel_size': 0.04, 'truncation': 6}
psf_class = PSF(**kwargs_psf)
lens_light_model_list = ['SERSIC_ELLIPSE']
lightModel = LightModel(lens_light_model_list)
supersampled_indexes = np.zeros((numPix, numPix), dtype=bool)
kwargs_numerics = {'supersampling_factor': 4}
imageModel = ImageModel(data_class, psf_class, lens_light_model_class=lightModel, kwargs_numerics=kwargs_numerics)

flux_mis_list = []
flux_truth_list = []
# for folder in folder_list:
for folder in [folder_list[0]]:    
    key = folder    
    kwargs_truth_i = result_dic[key][0]
    kwargs_truth_light = kwargs_truth_i['kwargs_lens_light']
    kwargs_result_i = result_dic[key][1]
    kwargs_result_light = kwargs_result_i['kwargs_lens_light']
    image_truth = imageModel.image(kwargs_lens_light=kwargs_truth_light, unconvolved=True)
    image_result = imageModel.image(kwargs_lens_light=kwargs_result_light, unconvolved=True)
    
    kwargs_truth_ps = kwargs_truth_i['kwargs_ps']
    x = kwargs_truth_ps['ra_image'] / 0.04 + int(numPix/2)
    y = kwargs_truth_ps['dec_image'] / 0.04 + int(numPix/2)
    plt.imshow(( abs(image_truth -image_result )/image_truth ), origin='lower')
    for i in range(len(x)):
        	plt.plot(x[i], y[i], 'xr')
        	# ax.text(x, x, 'x', fontsize=20, color='k') 
    flux_turth = np.array([image_truth[int(x[i]), int(y[i])] for i in range(len(x))])
    flux_result = np.array([image_result[int(x[i]), int(y[i])] for i in range(len(x))])
    flux_mis = flux_result - flux_turth
    flux_mis_list.append(flux_mis)
    flux_truth_list.append(flux_turth)
    plt.colorbar()
    plt.close()
#%%
mis = []
for i in range(len(folder_list)):
    for j in range(4):
        plt.scatter(i, flux_mis_list[i][j]/flux_truth_list[i][j])
        mis.append(flux_mis_list[i][j]/flux_truth_list[i][j])
plt.show()

plt.hist(mis)