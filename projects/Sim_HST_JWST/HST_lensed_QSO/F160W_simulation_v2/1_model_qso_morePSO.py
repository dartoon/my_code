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

folder_list = glob.glob('simulations_700_subg30/sim_lens_ID_subg30_7??')
folder_list.sort()

index = 0 #int(sys.argv[1])
print("Which index:", index)
# index = int(input("which index:\n"))
savename = 'PSFerr025_notPSFinter_morePSO_0.pkl'
# savename = 'PSFerr025_notPSFinter_2nd.pkl'
# savename = 'PSFerr025_notPSFinter.pkl'

#%%plot_spec_filter
save_pkl_folder = 'AGN_result_folder/'
    
for folder in [folder_list[index]]:
    ID = folder[-3:]
    folder = folder + '/'
    print(folder)
    model_lists, para_s, lens_info= pickle.load(open(folder+'sim_kwargs.pkl','rb'))
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
    save_file = save_pkl_folder+'idx{0}_ID'.format(index)+ID+'_'+savename
    if glob.glob(save_file) == []:    
        lens_data = pyfits.getdata(folder+'Drz_QSO_image.fits')
        len_std = pyfits.getdata(folder+'noise_map.fits')
        lens_mask = cr_mask(lens_data, 'normal_mask.reg')
        framesize = 81
        ct = int((len(lens_data) - framesize)/2)
        lens_data = lens_data[ct:-ct,ct:-ct]
        len_std = len_std[ct:-ct,ct:-ct]
        lens_mask = (1-lens_mask)[ct:-ct,ct:-ct]
        plt.imshow(lens_data * lens_mask, origin='lower',cmap='gist_heat', norm=LogNorm())
        plt.colorbar()
        exp_time = 599.* 2 * 8
        deltaPix = 0.08
        
        x, y =find_loc_max(lens_data)
        center = (len(lens_data)-1)/2
        QSO_pos= []
        for i in range(len(x)):
            QSO_pos.append([float(x[i]-0.5), float(y[i]-0.5)])
        xy_index = np.indices((len(lens_data),len(lens_data)))
        # for i in range(len(x)):

        count = 0
        for i in range(len(x)):
            x0, y0 =  (float(x[i]) - center) * deltaPix , (float(y[i]) - center) * deltaPix
#            print(x0, y0)
            ds = (x0- kwargs_ps['ra_image'] )**2 + (y0-kwargs_ps['dec_image'])**2
            if ds.min()<0.01:
#                print(np.where(ds == ds.min()), ds.min())
#                print("shift: ", kwargs_ps['ra_image'][ds == ds.min()] - x0, kwargs_ps['dec_image'][ds == ds.min()]-y0)
                kwargs_ps['ra_image'][ds == ds.min()] = x0
                kwargs_ps['dec_image'][ds == ds.min()] = y0
                count += 1
                
        if count!= len(kwargs_ps['ra_image']):
            raise ValueError("the PS positions is not assigned correctly")
        for i in range(len(kwargs_ps['ra_image'])):
            plt.plot(kwargs_ps['ra_image']/deltaPix+center, kwargs_ps['dec_image']/deltaPix+center, 'go')
        plt.close()         
        
        psf = pyfits.getdata(folder+'Drz_PSF.fits')
        psf_fsize = 77
        psf_half_r = int(psf_fsize/2)
        psf_peak = np.where(psf==psf.max())
        psf_peak = [psf_peak[0][0], psf_peak[1][0]]
        psf = psf[psf_peak[0]-psf_half_r:psf_peak[0]+psf_half_r+1,psf_peak[1]-psf_half_r:psf_peak[1]+psf_half_r+1]

        kwargs_data = sim_util.data_configure_simple(numPix = framesize, deltaPix = deltaPix) #,inverse=True)
        kwargs_data['image_data'] = lens_data
        kwargs_data['noise_map'] = len_std
        
        data_class = ImageData(**kwargs_data)
        psf_err_map = np.ones_like(psf)*0.25
        # psf_err_map[psf<0.003] = 0
        kwargs_psf = {'psf_type': 'PIXEL', 'kernel_point_source': psf, 'pixel_size': deltaPix, 'psf_error_map': psf_err_map}
        psf_class = PSF(**kwargs_psf)
        
        #%%
        # lens model choicers
        fixed_lens = []
        kwargs_lens_init = []
        kwargs_lens_sigma = []
        kwargs_lower_lens = []
        kwargs_upper_lens = []
        fixed_lens.append({}) 
        fixed_lens.append({'ra_0': 0, 'dec_0': 0})
        kwargs_lens_init = kwargs_lens_list
        kwargs_lens_sigma.append({'theta_E': .2, 'e1': 0.1, 'e2': 0.1, 'gamma': 0.1, 'center_x': 0.01, 'center_y': 0.01})
        kwargs_lower_lens.append({'theta_E': 0.01, 'e1': -0.5, 'e2': -0.5, 'gamma': kwargs_lens_init[0]['gamma']-0.5, 'center_x': -10, 'center_y': -10})
        kwargs_upper_lens.append({'theta_E': 10, 'e1': 0.5, 'e2': 0.5, 'gamma': kwargs_lens_init[0]['gamma']+0.5, 'center_x': 10, 'center_y': 10})
        kwargs_lens_sigma.append({'gamma1': 0.1, 'gamma2': 0.1})
        kwargs_lower_lens.append({'gamma1': -0.2, 'gamma2': -0.1})
        kwargs_upper_lens.append({'gamma1': 0.2, 'gamma2': 0.2})
        lens_params = [kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
        
        # lens light model choices
        fixed_lens_light = []
        kwargs_lens_light_init = []
        kwargs_lens_light_sigma = []
        kwargs_lower_lens_light = []
        kwargs_upper_lens_light = []
        kwargs_lens_light_init = kwargs_lens_light_list
        fixed_lens_light.append({})
        kwargs_lens_light_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.1, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
        kwargs_lower_lens_light.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.01, 'n_sersic': 0.5, 'center_x': -10, 'center_y': -10})
        kwargs_upper_lens_light.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10, 'n_sersic': 8, 'center_x': 10, 'center_y': 10})
        lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma, fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]
        
        # source model choices
        fixed_source = []
        kwargs_source_init = []
        kwargs_source_sigma = []
        kwargs_lower_source = []
        kwargs_upper_source = []
        fixed_source.append({})
        kwargs_source_init = kwargs_source_list
        kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.05, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.001, 'n_sersic': .5, 'center_x': -10, 'center_y': -10})
        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10, 'n_sersic': 5., 'center_x': 10, 'center_y': 10})
        source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
        
        # point source choices
        fixed_ps = [{}]
        kwargs_ps_init = [kwargs_ps]
#        kwargs_ps_init = kwargs_result['kwargs_ps']
        kwargs_ps_sigma = [{'ra_image': 0.01 * np.ones(len(kwargs_ps_init[0]['ra_image'])), 'dec_image': 0.01 * np.ones(len(kwargs_ps_init[0]['ra_image']))}]
        kwargs_lower_ps = [{'ra_image': -10 * np.ones(len(kwargs_ps_init[0]['ra_image'])), 'dec_image': -10 * np.ones(len(kwargs_ps_init[0]['ra_image']))}]
        kwargs_upper_ps = [{'ra_image': 10* np.ones(len(kwargs_ps_init[0]['ra_image'])), 'dec_image': 10 * np.ones(len(kwargs_ps_init[0]['ra_image']))}]
        ps_params = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]
        
        # Set cosmo
        fixed_special = {}
        kwargs_special_init = {}
        kwargs_special_sigma = {}
        kwargs_lower_special = {}
        kwargs_upper_special = {}
        kwargs_special_init = {'D_dt': TD_distance}
        kwargs_special_sigma = {'D_dt': 500}
        kwargs_lower_special = {'D_dt': TD_distance/2}
        kwargs_upper_special = {'D_dt': TD_distance*1.5}
        #Suggested by Simon: 
        astrometry_sigma = 0.005
        ximg, yimg = kwargs_ps_init[0]['ra_image'], kwargs_ps_init[0]['dec_image']
        kwargs_special_init['delta_x_image'], kwargs_special_init['delta_y_image'] = np.zeros_like(ximg), np.zeros_like(yimg)
        kwargs_special_sigma['delta_x_image'], kwargs_special_sigma['delta_y_image'] = np.ones_like(ximg) * astrometry_sigma, np.ones_like(yimg) * astrometry_sigma
        kwargs_lower_special['delta_x_image'], kwargs_lower_special['delta_y_image'] = np.ones_like(ximg) * (-1), np.ones_like(yimg) * (-1)
        kwargs_upper_special['delta_x_image'], kwargs_upper_special['delta_y_image'] = np.ones_like(ximg) * (1), np.ones_like(yimg) * (1)
        special_params = [kwargs_special_init, kwargs_special_sigma, fixed_special, kwargs_lower_special, kwargs_upper_special]
        
        kwargs_params = {'lens_model': lens_params,
                        'source_model': source_params,
                        'lens_light_model': lens_light_params,
                        'point_source_model': ps_params,
                        'special': special_params,
                        'point_source_offset': True}
        
        # numerical options and fitting sequences
        num_source_model = len(source_model_list)
        
        kwargs_likelihood = {'check_bounds': True,
                              'force_no_add_image': False,
                              'source_marg': False,
                              'image_position_uncertainty': astrometry_sigma,
                              'astrometric_likelihood': True,
                              'check_matched_source_position': True,
                              'source_position_tolerance': 0.01,
                              'time_delay_likelihood': True,
                              'image_likelihood_mask_list': [lens_mask]
                                      }
        kwargs_numerics = {'supersampling_factor': 3}
        image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
        multi_band_list = [image_band]
        kwargs_data_joint = {'multi_band_list': multi_band_list, 'multi_band_type': 'multi-linear',
                            'time_delays_measured': TD_obs[1:],
                            'time_delays_uncertainties': TD_err_l[1:]}
        
        kwargs_model = {'lens_model_list': lens_model_list, 
                          'lens_light_model_list': lens_light_model_list,
                          'source_light_model_list': source_model_list,
                        'point_source_model_list': point_source_list
                          }
        fitting_seq = FittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints,
                                      kwargs_likelihood, kwargs_params)
        
        start_time = time.time()
        # fitting_kwargs_list_0 = [
        #                         ['PSO', {'sigma_scale': 1., 'n_particles': 150, 'n_iterations': 200}],
        #                         ['PSO', {'sigma_scale': 1., 'n_particles': 200, 'n_iterations': 400}],
        #                         ['PSO', {'sigma_scale': 1., 'n_particles': 200, 'n_iterations': 400}]
        #                         ]
        
        # fitting_seq.fit_sequence(fitting_kwargs_list_0)
        
        kwargs_psf_iter = {'num_iter': 150, 'psf_iter_factor': 0.5,
                            'stacking_method': 'median', 
                            'keep_psf_error_map': True, 
                            'psf_symmetry': 1, 
                            'block_center_neighbour': 0.05}
        
        # # kwargs_psf_iter = {'num_iter': 150, 'psf_iter_factor': 0.5,
        # #                     'stacking_method': 'median', 
        # #                     'keep_psf_error_map': False, 
        # #                     'psf_symmetry': 1, 
        # #                     'block_center_neighbour': 0.05}        
        
        fitting_kwargs_list_1 = [
                        ['PSO', {'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}],
                        # ['psf_iteration', kwargs_psf_iter],
                        ['PSO', {'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}],
                        # ['psf_iteration', kwargs_psf_iter],
                        ['PSO', {'sigma_scale': 1., 'n_particles': 150, 'n_iterations': 100}],
                        # ['psf_iteration', kwargs_psf_iter],
                        ['PSO', {'sigma_scale': .1, 'n_particles': 150, 'n_iterations': 100}],
                        # ['psf_iteration', kwargs_psf_iter],
                        ['PSO', {'sigma_scale': 1., 'n_particles': 200, 'n_iterations': 200}],
                        # ['psf_iteration', kwargs_psf_iter],
                        ['PSO', {'sigma_scale': .1, 'n_particles': 200, 'n_iterations': 200}],
                        # ['psf_iteration', kwargs_psf_iter],
                        ['MCMC', {'n_burn': 500, 'n_run': 600, 'walkerRatio': 6, 'sigma_scale': 0.1}],
                        ]
        chain_list = fitting_seq.fit_sequence(fitting_kwargs_list_1)
        kwargs_result_best = fitting_seq.best_fit()
        end_time = time.time()
        print(end_time - start_time, 'total time needed for computation')
        print('============ CONGRATULATION, YOUR JOB WAS SUCCESSFUL ================ ')
        #Save in pickle
        fix_setting =  [fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_special]
#        cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.)  #!!!Wrong cosmos
#        td_cosmo = TDCosmography(z_l, z_s, kwargs_model, cosmo_fiducial=cosmo)
        # make instance of parameter class with given model options, constraints and fixed parameters #
        param = Param(kwargs_model, fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_special, 
                      kwargs_lens_init=kwargs_result_best['kwargs_lens'], **kwargs_constraints)
        sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[-1]
        mcmc_new_list = []
        labels_new = [r"$\gamma$", r"$D_{\Delta t}$","H$_0$" ]
        for i in range(len(samples_mcmc)):
            # transform the parameter position of the MCMC chain in a lenstronomy convention with keyword arguments #
            kwargs_result = param.args2kwargs(samples_mcmc[i])
            D_dt = kwargs_result['kwargs_special']['D_dt']
#            fermat_pot = td_cosmo.fermat_potential(kwargs_result['kwargs_lens'], kwargs_result['kwargs_ps'])
        #    delta_fermat_12 = fermat_pot[0] - fermat_pot[2]
            gamma = kwargs_result['kwargs_lens'][0]['gamma']
        #    phi_ext, gamma_ext = kwargs_result['kwargs_lens'][1]['gamma1'], kwargs_result['kwargs_lens'][1]['gamma2']
            mcmc_new_list.append([gamma, D_dt, cal_h0(z_l ,z_s, D_dt)])        
        
        pickle.dump([multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list], open(save_file, 'wb'))
    #%%Print fitting result:
    multi_band_list, kwargs_model, kwargs_result_best, chain_list, fix_setting, mcmc_new_list = pickle.load(open(save_file,'rb'))
    # fixed_lens, fixed_source, fixed_lens_light, fixed_ps, fixed_cosmo = fix_setting
    labels_new = [r"$\gamma$", r"$D_{\Delta t}$","H$_0$" ]
    modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result_best, arrow_size=0.02, cmap_string="gist_heat")
    f, axes = modelPlot.plot_main()
    f.show()
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
    plt.show()
