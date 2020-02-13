#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 11:02:31 2019

@author: Dartoon

Fit the Seed 211 - 216 with a simple model
"""
# import the necessary python modules
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

steps = [1000, 10000]
H0_true = 70.65595
for fix_gamma in [False]:
    for seed in [211]:
        print("Runing seed:", seed, "fix_gamma:", fix_gamma, "MC steps:", steps)
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
        kwargs_lens_sigma.append({'theta_E': .1, 'e1': 0.1, 'e2': 0.1, 'gamma': 0.1, 'center_x': 0.005, 'center_y': 0.005}) #previous: 'center_x': 0.1, 'center_y': 0.1
        kwargs_lower_lens.append({'theta_E': 0.01, 'e1': -0.5, 'e2': -0.5, 'gamma': 1.5, 'center_x': 0-0.01, 'center_y': 0-0.01})
        kwargs_upper_lens.append({'theta_E': 10, 'e1': 0.5, 'e2': 0.5, 'gamma': 2.5, 'center_x': 0+0.01, 'center_y': 0+0.01})
        
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
                              'solver_type': 'PROFILE',  # 'PROFILE_SHEAR', 'NONE', # any proposed lens model must satisfy the image positions appearing at the position of the point sources being sampeld
                              'Ddt_sampling': time_delay_likelihood,  # sampling of the time-delay distance                      
                             }
        
        # we can define un-correlated Gaussian priors on specific parameters explicitly
        # e.g. power-law mass slope of the main deflector
        #prior_lens = [[0, 'gamma', 2, 0.15]] # [[index_model, 'param_name', mean, 1-sigma error], [...], ...]
        # e.g. source size of the emission region
        #prior_special = []
            
        kwargs_likelihood = {  
                             'image_position_uncertainty': astrometry_sigma,  # astrometric uncertainty of image positions
                             'source_position_likelihood': True,  # evaluates how close the different image positions match the source positons
                             'image_position_likelihood': True, # evaluate point source likelihood given the measured image positions
                             'time_delay_likelihood': time_delay_likelihood,  # evaluating the time-delay likelihood
                             #'prior_lens': prior_lens,
        #                     'prior_special': prior_special,
                             'check_matched_source_position': True, # check non-linear solver and disgard non-solutions
                             'source_position_tolerance': 0.01,
                             'check_bounds': True,  # check parameter bounds and punish them
                             'source_position_sigma': 0.001
                            }
        
        from lenstronomy.Workflow.fitting_sequence import FittingSequence
        fitting_seq = FittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
        
        fitting_kwargs_list = [#['update_settings', {'lens_add_fixed': [[0, ['gamma']]]}],  # you can add additional fixed parameters if you want
                               ['PSO', {'sigma_scale': 1., 'n_particles': 300, 'n_iterations': 300}],
                               ['PSO', {'sigma_scale': 1., 'n_particles': 300, 'n_iterations': 200}],
                               ['PSO', {'sigma_scale': 1., 'n_particles': 300, 'n_iterations': 200}]
                               ]
        
        start_time = time.time()
        chain_list_pso = fitting_seq.fit_sequence(fitting_kwargs_list)
        kwargs_result = fitting_seq.best_fit()
        end_time = time.time()
        print(end_time - start_time, 'total time needed for computation')
        print('============ CONGRATULATION, YOUR JOB WAS SUCCESSFUL ================ ')
        
        #%%
        kwargs_result = fitting_seq.best_fit(bijective=True)
    
        from lenstronomy.Plots import chain_plot as chain_plot
        for i in range(len(chain_list_pso)):
            chain_plot.plot_chain_list(chain_list_pso, i)
            plt.close()
            
        #and now we run the MCMC
        fitting_kwargs_list = [
            ['MCMC', {'n_burn': steps[0], 'n_run': steps[1], 'walkerRatio': 30,'sigma_scale': 0.1}]
        ]
        chain_list_mcmc = fitting_seq.fit_sequence(fitting_kwargs_list)
        kwargs_result = fitting_seq.best_fit()    
    
        args_result = fitting_seq.param_class.kwargs2args(**kwargs_result)
        logL = fitting_seq.likelihoodModule.logL(args_result, verbose=True)
        
        sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list_mcmc[0]
        
        # import the parameter handling class #
        from lenstronomy.Sampling.parameters import Param
        import lenstronomy.Util.param_util as param_util
        param = Param(kwargs_model, fixed_lens, kwargs_fixed_ps=fixed_ps, kwargs_fixed_special=fixed_special, 
                      kwargs_lens_init=kwargs_result['kwargs_lens'], **kwargs_constraints)
        # the number of non-linear parameters and their names #
        num_param, param_list = param.num_param()
        
        lensModel = LensModel(kwargs_model['lens_model_list'])
        lensModelExtensions = LensModelExtensions(lensModel=lensModel) 
        
        mcmc_new_list = []
        if fix_gamma:
            labels_new = [r"$\theta_E$", r"$\phi_{lens}$", r"$q$",r"$D_{dt}$", r"$H_0$"]
        else:
            labels_new = [r"$\theta_E$", r"$\gamma$", r"$\phi_{lens}$", r"$q$",r"$D_{dt}$", r"$H_0$"]
            
    #        for i in range(len(samples_mcmc)):
        trans_steps = int(np.min([len(samples_mcmc)/10, 40000]))
        for i in range(trans_steps):
            # transform the parameter position of the MCMC chain in a lenstronomy convention with keyword arguments #
            kwargs_out = param.args2kwargs(samples_mcmc[-trans_steps+i])
            kwargs_lens_out, kwargs_special_out, kwargs_ps_out = kwargs_out['kwargs_lens'], kwargs_out['kwargs_special'], kwargs_out['kwargs_ps']
            
            # compute 'real' image position adding potential astrometric shifts
            x_pos, y_pos = kwargs_ps_out[0]['ra_image'], kwargs_ps_out[0]['dec_image']
            
            # extract quantities of the main deflector
            theta_E = kwargs_lens_out[0]['theta_E']
            e1, e2 = kwargs_lens_out[0]['e1'], kwargs_lens_out[0]['e2']
            phi, q = param_util.ellipticity2phi_q(e1, e2)
            if fix_gamma:
                new_chain = [theta_E, phi, q]
            else:
                gamma = kwargs_lens_out[0]['gamma']
                new_chain = [theta_E, gamma, phi, q]
            D_dt = kwargs_special_out['D_dt']
            new_chain.append(D_dt)
            new_chain.append(cal_h0(z_lens ,z_source, D_dt))
            #source_size = np.random.uniform(high=1, low=0)
            mcmc_new_list.append(np.array(new_chain))
            if i/2000 > (i-1)/2000 :
                print("total",len(samples_mcmc), "finished translate:", i)
        
        # plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True,
        #                      quantiles=[0.16, 0.5, 0.84],
        #                      title_kwargs={"fontsize": 15}, label_kwargs = {"fontsize": 25},
        #                      levels=1.0 - np.exp(-0.5 * np.array([1.,2.]) ** 2))
        
        import os
        if os.path.exists('./MCsteps{0}_{1}'.format(steps[0], steps[1]))==False:
            os.mkdir('./MCsteps{0}_{1}'.format(steps[0], steps[1]))
        
        if fix_gamma:
            # plot.savefig("MCsteps{0}_{1}/corner_plot_SIE#{2}.pdf".format(steps[0], steps[1], seed-210))
            datafile = 'MCsteps{0}_{1}/run_result_SIE.txt'.format(steps[0], steps[1])
            picklename = 'MCsteps{0}_{1}/sampler_results_SIE#{2}.pkl'.format(steps[0], steps[1], seed-210)
        else:
            # plot.savefig("MCsteps{0}_{1}/corner_plot_SPEMD#{2}.pdf".format(steps[0], steps[1], seed-210))                     
            datafile = 'MCsteps{0}_{1}/run_result_SPEMD.txt'.format(steps[0], steps[1])
            picklename = 'MCsteps{0}_{1}/sampler_results_SPEMD#{2}.pkl'.format(steps[0], steps[1], seed-210)
        plt.close()
        
        mcmc_new_list = np.asarray(mcmc_new_list)
        if_file = glob.glob(datafile)
        if if_file == []:
            para_result =  open(datafile,'w') 
        else:
            para_result = open(datafile,'r+')
            para_result.read()
        para_result.write("seed = {0}: ".format(seed-210))    
        para_result.write('H_0: {0:.3f}pm{1:.3f}, total_Chisq:{2:.3f}, '.format(np.median(mcmc_new_list[:,-1]), 
                          (np.percentile(mcmc_new_list[:,-1],84,axis=0) - np.percentile(mcmc_new_list[:,-1],16,axis=0))/2.
                          , -logL* 2))
        if not fix_gamma:
            para_result.write('gamma: {0:.3f}pm{1:.3f}, '.format(np.median(mcmc_new_list[:,1]),
                              (np.percentile(mcmc_new_list[:,1],84,axis=0) - np.percentile(mcmc_new_list[:,1],16,axis=0))/2.))
        para_result.write('q: {0:.3f}pm{1:.3f} \n'.format(np.median(mcmc_new_list[:,-3]), 
                          (np.percentile(mcmc_new_list[:,-3],84,axis=0) - np.percentile(mcmc_new_list[:,-3],16,axis=0))/2.))
        para_result.close()    
        if len(chain_list_mcmc[0][1])>1200000:
            chain_list_mcmc[0][1] = chain_list_mcmc[0][1][-1200000:]
        fit_result = [kwargs_result, chain_list_mcmc, chain_list_pso]
        trans_result = [mcmc_new_list, labels_new]
        pickle.dump([fit_result, trans_result], open(picklename, 'wb'))            
        
