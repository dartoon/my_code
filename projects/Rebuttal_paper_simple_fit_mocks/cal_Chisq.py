#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 14:51:56 2019

@author: Dartoon

Calculate the Chisq
"""
import numpy as np
import time
import corner

from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.lens_model_extensions import LensModelExtensions
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from lenstronomy.Cosmo.lens_cosmo import LensCosmo
from lenstronomy.Util import constants
from astropy.cosmology import FlatLambdaCDM
from  lenstronomy.Plots import lens_plot
import glob
import matplotlib.pyplot as plt
import pickle

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

H0_true = 70.65595
#noise_seed = 2

#folder = 'MCsteps[1000, 10000]_first_run_linear_True/'
#folder = 'MCsteps1000_10000_second_run_linear_True/'
folder = 'MCsteps1000_10000_deflector_10mas/'

for fix_gamma in [False]:
    for seed in [216]:
        if seed == 211:
            z_lens ,z_source = [0.640, 2.408]
            ximg = [0.41530, 0.17330, 1.22509, -1.18836]  # image positions in relative RA (arc seconds)
            yimg = [-1.19026, 1.24663, 0.14337, -0.11013]  # image positions in relative DEC (arc seconds)
            d_dt_sigma = np.array([0.250, 0.250, 0.250])  # 1-sigma uncertainties in the time-delay measurement (in units of days)
            d_dt_measured = [0.25322, 4.15310, 8.87434]
            D_dt_true = 3898
            kwargs_lens_list = [{'theta_E': 1.23651, 'center_x': 0, 'center_y': 0, 'e1': 0.05207, 'gamma': 2.00000, 'e2': 0.01075}]
            e1, e2 = kwargs_lens_list[0]['e1'], kwargs_lens_list[0]['e2']
        elif seed == 212:
            z_lens ,z_source = [0.54587, 1.94048]
            ximg = [0.30349, 0.56195, 1.09881, -1.06851]	  # image positions in relative RA (arc seconds)
            yimg = [1.14822, -1.00578, -0.32878, -0.20075]  # image positions in relative DEC (arc seconds)
            d_dt_sigma = np.array([0.250, 0.250, 0.250])  # 1-sigma uncertainties in the time-delay measurement (in units of days)
            d_dt_measured = [3.63535, 4.21165, 10.68162]
            D_dt_true = 3367
            kwargs_lens_list = [{'theta_E': 1.14303, 'center_x': 0, 'center_y': 0, 'e1': 0.05884, 'gamma': 2.00000, 'e2': 0.00370}]
            e1, e2 = kwargs_lens_list[0]['e1'], kwargs_lens_list[0]['e2']
        elif seed == 213:
            z_lens ,z_source = [0.52697, 2.23751]
            ximg = [0.88708, -0.97571, -0.75385, 0.64415]	  # image positions in relative RA (arc seconds)
            yimg = [0.84356, -0.70413, 0.80061, -0.85652]  # image positions in relative DEC (arc seconds)
            d_dt_sigma = np.array([0.250, 0.250, 0.250])  # 1-sigma uncertainties in the time-delay measurement (in units of days)
            d_dt_measured = [1.52899, 5.92121, 6.86775]
            D_dt_true = 3059
            kwargs_lens_list =  [{'theta_E': 1.79689, 's_scale': 0.31301, 'center_x': 0, 'center_y': 0, 'e1': -0.01025, 'gamma': 2.00000, 'e2': -0.05716}]
            e1, e2 = kwargs_lens_list[0]['e1'], kwargs_lens_list[0]['e2']
        elif seed == 214:
            z_lens ,z_source = [0.47201, 1.82471]
            ximg = [1.12211, -0.72429, -0.03085, -0.37896]	  # image positions in relative RA (arc seconds)
            yimg = [0.01948, 0.76705, 1.02654, -0.79019]  # image positions in relative DEC (arc seconds)
            d_dt_sigma = np.array([0.250, 0.250, 0.250])  # 1-sigma uncertainties in the time-delay measurement (in units of days)
            d_dt_measured = [3.08772, 3.72799, 8.89355]
            D_dt_true = 2807
            kwargs_lens_list = [{'theta_E': 1.68307, 's_scale': 0.31427, 'center_x': 0, 'center_y': 0, 'e1': -0.04574, 'gamma': 2.00000, 'e2': 0.03103}]
            e1, e2 = kwargs_lens_list[0]['e1'], kwargs_lens_list[0]['e2']
        elif seed == 215:
            z_lens ,z_source = [0.49269, 1.51525]
            ximg = [-0.53603, 0.74086, 0.84869, -0.70392]	  # image positions in relative RA (arc seconds)
            yimg = [0.96151, -0.76043, 0.49576, -0.59825]  # image positions in relative DEC (arc seconds)
            d_dt_sigma = np.array([0.250, 0.250, 0.271])  # 1-sigma uncertainties in the time-delay measurement (in units of days)
            d_dt_measured = [3.42016, 9.03260, 13.55273]
            D_dt_true = 3213
            kwargs_lens_list = [{'Rs': 1.48138, 'center_x': 0.00000, 'center_y': 0.00000, 'e1': 0.01711, 'sigma0': 0.80000, 'e2': 0.04999}, 
                            {'alpha_Rs': 1.00000, 'Rs': 31.18490, 'center_x': 0.00000, 'center_y': 0.00000, 'e1': 0.01345, 'e2': 0.05118}]
            e1, e2 = kwargs_lens_list[0]['e1'], kwargs_lens_list[0]['e2']
        elif seed == 216:
            z_lens ,z_source = [0.29370, 2.26521]
            ximg = [0.72658, 0.53087, 0.96278, -0.73249]	  # image positions in relative RA (arc seconds)
            yimg = [1.26847, -1.13426, -0.69393, 0.00453]  # image positions in relative DEC (arc seconds)
            d_dt_sigma = np.array([0.250, 0.250, 0.250])  # 1-sigma uncertainties in the time-delay measurement (in units of days)
            d_dt_measured = [5.26037, 4.92182, 12.11084]
            D_dt_true = 1472
            kwargs_lens_list = [{'Rs': 1.96876, 'center_x': 0.00000, 'center_y': 0.00000, 'e1': 0.05257, 'sigma0': 0.20000, 'e2': -0.02539}, 
                                {'alpha_Rs': 4.50000, 'Rs': 34.49701, 'center_x': 0.00000, 'center_y': 0.00000, 'e1': 0.04156, 'e2': -0.01512}]
            e1, e2 = kwargs_lens_list[0]['e1'], kwargs_lens_list[0]['e2']
            
        astrometry_sigma = 0.004  # 1-sigma astrometric uncertainties of the image positions (assuming equal precision for all images in RA/DEC directions)
        #np.random.seed(noise_seed)
        ximg_measured = ximg #+ np.random.normal(0, astrometry_sigma, 4)
        yimg_measured = yimg #+ np.random.normal(0, astrometry_sigma, 4)
        
        # here we create a keyword list with all the data elements. If you only have partial information about your lens,
        # only provide the quantities you have.
        kwargs_data_joint = {'time_delays_measured': d_dt_measured,
                             'time_delays_uncertainties': d_dt_sigma,
                             'ra_image_list': [ximg_measured], 'dec_image_list': [yimg_measured]}
        
        # ==================
        # lens model choices
        # ==================
        lens_model_list = ['SPEMD']
        fixed_lens = []
        kwargs_lens_init = []
        kwargs_lens_sigma = []
        kwargs_lower_lens = []
        kwargs_upper_lens = []
        # SPEMD parameters
        if fix_gamma:
            fixed_lens.append({'gamma': 2.0}) 
        else:
            fixed_lens.append({})
            
        # initial parameter guess
        #kwargs_lens_init.append(kwargs_lens[0])
        kwargs_lens_init.append({'theta_E': 1.20, 'gamma': 2., 'center_x': 0, 'center_y': 0, 'e1': e1, 'e2': e2})
        kwargs_lens_sigma.append({'theta_E': .1, 'e1': 0.1, 'e2': 0.1, 'gamma': 0.1, 'center_x': 0.1, 'center_y': 0.1})
        kwargs_lower_lens.append({'theta_E': 0.01, 'e1': -0.5, 'e2': -0.5, 'gamma': 1., 'center_x': e1-0.08, 'center_y': e2-0.08})
        kwargs_upper_lens.append({'theta_E': 10, 'e1': 0.5, 'e2': 0.5, 'gamma': 3, 'center_x': e1+0.08, 'center_y': e2+0.08})
        
        # combine all parameter options for lenstronomy
        lens_params = [kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
        
        # =========================
        # image position parameters
        # =========================
        # we chose to model the image positions in the lensed plane (we know where they appear)
        point_source_list = ['LENSED_POSITION']
        # We fix the image position coordinates.
        fixed_ps = [{}]  # we fix the image position coordinates
        # these lines below actually don't matter when you keep the image position fixed
        kwargs_ps_init = [{'ra_image': ximg_measured, 'dec_image': yimg_measured}] # the initial guess for the appearing image positions is: at the image position.
        kwargs_ps_sigma = [{'ra_image': 0.01 * np.ones(len(ximg)), 'dec_image': 0.01 * np.ones(len(ximg))}]
        kwargs_lower_ps = [{'ra_image': -10 * np.ones(len(ximg)), 'dec_image': -10 * np.ones(len(ximg))}]
        kwargs_upper_ps = [{'ra_image': 10* np.ones(len(ximg)), 'dec_image': 10 * np.ones(len(ximg))}]
        
        # combine all parameter options for lenstronomy
        ps_params = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]
        
        fixed_special = {}
        kwargs_special_init = {}
        kwargs_special_sigma = {}
        kwargs_lower_special = {}
        kwargs_upper_special = {}
        # ===================
        # Time-delay distance
        # ===================
        # with time-delay information, we can measure the time-delay distance (units physical Mpc)
        # if you want to fix the cosmology and instead use the time-delay information to constrain the lens model, out-comment the line below
        kwargs_special_init['D_dt'] = D_dt_true
        kwargs_special_sigma['D_dt'] = 1000
    #    kwargs_lower_special['D_dt'] = D_dt_true*H0_true/100.
    #    kwargs_upper_special['D_dt'] = D_dt_true*H0_true/50.
        kwargs_lower_special['D_dt'] = 0
        kwargs_upper_special['D_dt'] = 20000    
        
        special_params = [kwargs_special_init, kwargs_special_sigma, fixed_special, kwargs_lower_special, kwargs_upper_special]
        
        # combined parameter settings
        kwargs_params = {'lens_model': lens_params,
                        'point_source_model': ps_params,
                        'special': special_params}
        
        # our model choices
        kwargs_model = {'lens_model_list': lens_model_list, 
                        'point_source_model_list': point_source_list
                         }
        
        time_delay_likelihood = True  # bool, set this True or False depending on whether time-delay information is available and you want to make use of its information content.
        image_position_likelihood = True  # bool, evaluating the image position likelihood (in combination with astrometric errors)
        
        kwargs_constraints = {'num_point_source_list': [len(ximg)],  
                              'solver_type': 'NONE',  # 'PROFILE_SHEAR', 'NONE', # any proposed lens model must satisfy the image positions appearing at the position of the point sources being sampeld
                              'Ddt_sampling': time_delay_likelihood,  # sampling of the time-delay distance                      
                             }
        
        # we can define un-correlated Gaussian priors on specific parameters explicitly
        # e.g. power-law mass slope of the main deflector
        #prior_lens = [[0, 'gamma', 2, 0.15]] # [[index_model, 'param_name', mean, 1-sigma error], [...], ...]
        # e.g. source size of the emission region
        #prior_special = []
            
        kwargs_likelihood = {  
                             'position_uncertainty': astrometry_sigma,  # astrometric uncertainty of image positions
                             'source_position_likelihood': True,  # evaluates how close the different image positions match the source positons
                             'image_position_likelihood': True, # evaluate point source likelihood given the measured image positions
                             'time_delay_likelihood': time_delay_likelihood,  # evaluating the time-delay likelihood
                             #'prior_lens': prior_lens,
        #                     'prior_special': prior_special,
                             'check_solver': True, # check non-linear solver and disgard non-solutions
                             'solver_tolerance': 0.001,
                             'check_bounds': True,  # check parameter bounds and punish them
                            }
        
        from lenstronomy.Workflow.fitting_sequence import FittingSequence
        fitting_seq = FittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
        if fix_gamma == True:
            result = pickle.load(open(folder+'sampler_results_SIE#{0}.pkl'.format(seed-210),'rb'))
        elif fix_gamma == False:
            result = pickle.load(open(folder+'sampler_results_SPEMD#{0}.pkl'.format(seed-210),'rb'))
        fit_result, trans_result = result 
        kwargs_result, chain_list_mcmc, chain_list_pso = fit_result
        args_result = fitting_seq.param_class.kwargs2args(**kwargs_result)
        mcmc_new_list, labels_new = trans_result
        logL = fitting_seq.likelihoodModule.logL(args_result, verbose=True)
#        print "note:", seed-210, fix_gamma, round(-logL*2,3)
        if fix_gamma == True:
            print "seed = {0}: ".format(seed-210), 'H_0: {0:.3f}pm{1:.3f}, total_Chisq:{2:.3f}, '.format(np.median(mcmc_new_list[:,-1]), 
                          (np.percentile(mcmc_new_list[:,-1],84,axis=0) - np.percentile(mcmc_new_list[:,-1],16,axis=0))/2.
                          , -logL* 2),  'q: {0:.3f}pm{1:.3f} \n'.format(np.median(mcmc_new_list[:,-3]), 
                          (np.percentile(mcmc_new_list[:,-3],84,axis=0) - np.percentile(mcmc_new_list[:,-3],16,axis=0))/2.)
        elif fix_gamma == False:
            print "seed = {0}: ".format(seed-210), 'H_0: {0:.3f}pm{1:.3f}, total_Chisq:{2:.3f}, '.format(np.median(mcmc_new_list[:,-1]), 
                          (np.percentile(mcmc_new_list[:,-1],84,axis=0) - np.percentile(mcmc_new_list[:,-1],16,axis=0))/2.
                          , -logL* 2), 'gamma: {0:.3f}pm{1:.3f}, '.format(np.median(mcmc_new_list[:,1]),
                              (np.percentile(mcmc_new_list[:,1],84,axis=0) - np.percentile(mcmc_new_list[:,1],16,axis=0))/2.),  'q: {0:.3f}pm{1:.3f} \n'.format(np.median(mcmc_new_list[:,-3]), 
                          (np.percentile(mcmc_new_list[:,-3],84,axis=0) - np.percentile(mcmc_new_list[:,-3],16,axis=0))/2.)            