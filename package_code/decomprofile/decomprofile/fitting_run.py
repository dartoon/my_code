#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 09:06:14 2020

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import time
import corner
import pickle
import copy


class FittingRun(object):
    """
    A class to perform the fitting task:
        - define the way to fitting: PSO and MCMC
        - save all the useful fitting materials, if assign save_pkl
    """
    def __init__(self, fitting_specify_class):
        self.fitting_specify_class = fitting_specify_class
        self.fitting_seq = fitting_specify_class.fitting_seq
    
    def fitting_kwargs(self, algorithm_list = ['PSO', 'MCMC'], setting_list = [None, None]):
        if len(algorithm_list) != len(setting_list):
            raise ValueError("The algorithm_list and setting_list should be in the same length.") 
        fitting_kwargs_list = []
        for i in range(len(algorithm_list)):
            if setting_list[i] is None:
                setting = fitting_setting_temp(algorithm_list[i])
            else:
                setting = setting_list[i]
            fitting_kwargs_list.append([algorithm_list[i], setting])
        self.fitting_kwargs_list = fitting_kwargs_list
    
    def run(self):
        self.fitting_kwargs(algorithm_list = ['PSO', 'MCMC'], setting_list = [None, None])
        fitting_specify_class = self.fitting_specify_class
        start_time = time.time()
        chain_list = self.fitting_seq.fit_sequence(self.fitting_kwargs_list)
        kwargs_result = self.fitting_seq.best_fit()
        ps_result = kwargs_result['kwargs_ps']
        source_result = kwargs_result['kwargs_source']
        if self.fitting_kwargs_list[-1][0] == 'MCMC':
            self.sampler_type, self.samples_mcmc, self.param_mcmc, self.dist_mcmc  = chain_list[-1]    
        end_time = time.time()
        print(end_time - start_time, 'total time spent for this computation (s)')
        print('============ CONGRATULATION, YOUR JOB WAS SUCCESSFUL ================ ')
        
        from lenstronomy.ImSim.image_linear_solve import ImageLinearFit
        imageLinearFit = ImageLinearFit(data_class=fitting_specify_class.data_class, psf_class=fitting_specify_class.psf_class,
                                        source_model_class=fitting_specify_class.lightModel,
                                        point_source_class=fitting_specify_class.pointSource, 
                                        kwargs_numerics=fitting_specify_class.kwargs_numerics)    
        image_reconstructed, error_map, _, _ = imageLinearFit.image_linear_solve(kwargs_source=source_result,
                                                                                 kwargs_ps=ps_result)
        from lenstronomy.Plots.model_plot import ModelPlot
        # this is the linear inversion. The kwargs will be updated afterwards
        modelPlot = ModelPlot(fitting_specify_class.kwargs_data_joint['multi_band_list'],
                              fitting_specify_class.kwargs_model, kwargs_result,
                              arrow_size=0.02, cmap_string="gist_heat", 
                              likelihood_mask_list=fitting_specify_class.kwargs_likelihood['image_likelihood_mask_list'] )    

        imageModel = fitting_specify_class.imageModel
        image_host_list = []  #The linear_solver before and after LensModelPlot could have different result for very faint sources.
        for i in range(len(source_result)):
            image_host_list.append(imageModel.source_surface_brightness(source_result, de_lensed=True,unconvolved=False,k=i))
        
        image_ps_list = []
        for i in range(len(ps_result)):
            image_ps_list.append(imageModel.point_source(ps_result, k = i))
        
        self.chain_list = chain_list
        self.kwargs_result = kwargs_result
        self.ps_result = ps_result
        self.source_result = source_result
        self.modelPlot = modelPlot
        self.imageLinearFit = imageLinearFit
        self.reduced_Chisq =  imageLinearFit.reduced_chi2(image_reconstructed, error_map)
        self.image_host_list = image_host_list
        self.image_ps_list = image_ps_list

    def run_diag(self, diag_list = None):
        from lenstronomy.Plots import chain_plot
        if diag_list is None:
            for i in range(len(self.chain_list)):
                f, axes = chain_plot.plot_chain_list(self.chain_list,i)
        else:
            for i in diag_list:
                f, axes = chain_plot.plot_chain_list(self.chain_list,i)
        plt.show()

    def model_plot(self):
        f, axes = plt.subplots(3, 3, figsize=(16, 16), sharex=False, sharey=False)
        self.modelPlot.data_plot(ax=axes[0,0], text="Data")
        self.modelPlot.model_plot(ax=axes[0,1])
        self.modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
        
        self.modelPlot.decomposition_plot(ax=axes[1,0], text='Host galaxy', source_add=True, unconvolved=True)
        self.modelPlot.decomposition_plot(ax=axes[1,1], text='Host galaxy convolved', source_add=True)
        self.modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
        
        self.modelPlot.subtract_from_data_plot(ax=axes[2,0], text='Data - Point Source', point_source_add=True)
        self.modelPlot.subtract_from_data_plot(ax=axes[2,1], text='Data - host galaxy', source_add=True)
        self.modelPlot.subtract_from_data_plot(ax=axes[2,2], text='Data - host galaxy - Point Source', source_add=True, point_source_add=True)
        f.tight_layout()
        plt.show()   
        
    def params_corner_plot(self):
        if self.fitting_kwargs_list[-1][0] == 'MCMC':
            samples_mcmc = self.samples_mcmc
            n, num_param = np.shape(samples_mcmc)
            plot = corner.corner(samples_mcmc, labels=self.param_mcmc, show_titles=True)
            plt.show()         
        else:
            print("Are you sure you have perform MCMC in the last fitting?")

    def flux_corner_plot(self):
        from lenstronomy.Sampling.parameters import Param
        fitting_specify_class = self.fitting_specify_class
        if self.fitting_kwargs_list[-1][0] == 'MCMC':
            param = Param(fitting_specify_class.kwargs_model, kwargs_fixed_source=fitting_specify_class.source_params[2],
                          kwargs_fixed_ps=fitting_specify_class.ps_params[2], **fitting_specify_class.kwargs_constraints)
            mcmc_new_list = []
            if len(fitting_specify_class.point_source_list) >0 :
                qso_labels_new = ["Quasar_{0} flux".format(i) for i in range(len(fitting_specify_class.point_source_list))]
                galaxy_labels_new = ["Galaxy_{0} flux".format(i) for i in range(len(fitting_specify_class.light_model_list))]
                labels_new = qso_labels_new + galaxy_labels_new
            else:
                labels_new = ["Galaxy_{0} flux".format(i) for i in range(len(fitting_specify_class.light_model_list))]
            if len(self.samples_mcmc) > 10000:  #Only save maximum 10000 chain results.
                trans_steps = [len(self.samples_mcmc)-10000, len(self.samples_mcmc)]
            else:
                trans_steps = [0, len(self.samples_mcmc)]
            for i in range(trans_steps[0], trans_steps[1]):
                kwargs_out = param.args2kwargs(self.samples_mcmc[i])
                kwargs_light_source_out = kwargs_out['kwargs_source']
                kwargs_ps_out =  kwargs_out['kwargs_ps']
                image_reconstructed, _, _, _ = self.imageLinearFit.image_linear_solve(kwargs_source=kwargs_light_source_out,
                                                                                      kwargs_ps=kwargs_ps_out)
                flux_list_quasar = []
                if len(fitting_specify_class.point_source_list) > 0:
                    for j in range(len(fitting_specify_class.point_source_list)):
                        image_ps_j = fitting_specify_class.imageModel.point_source(kwargs_ps_out, k=j)
                        flux_list_quasar.append(np.sum(image_ps_j))
                flux_list_galaxy = []
                for j in range(len(fitting_specify_class.light_model_list)):
                    image_j = fitting_specify_class.imageModel.source_surface_brightness(kwargs_light_source_out,unconvolved= False, k=j)
                    flux_list_galaxy.append(np.sum(image_j))
                mcmc_new_list.append(flux_list_quasar + flux_list_galaxy )
                if int(i/1000) > int((i-1)/1000) :
                    print(trans_steps[1]-trans_steps[0],
                          "MCMC samplers in total, finished translate:", i-trans_steps[0] )
            plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
            plt.show()        
            self.mcmc_new_list = mcmc_new_list
            self.labels_new = labels_new

    # def dump_result(self, tag = 'result'):
    #     fitting_specify_class = self.fitting_specify_class
    #     if hasattr(self, 'mcmc_new_list'):
    #         trans_paras = [self.mcmc_new_list, self.labels_new, 'mcmc_new_list, labels_new']
    #     else:
    #         trans_paras = []
    #     picklename= tag + '.pkl'
    #     best_fit = [self.source_result, self.image_host_list, self.ps_result, self.image_ps_list,
    #                 ['source_result', 'image_host_list', 'ps_result', 'image_ps_list'] ]
    #     chain_list_result = [self.chain_list, 'chain_list']
    #     kwargs_fixed_source=fitting_specify_class.source_params[2]
    #     kwargs_fixed_ps=fitting_specify_class.ps_params[2]
    #     model_classes = [fitting_specify_class.data_class, fitting_specify_class.psf_class,
    #                      fitting_specify_class.lightModel, fitting_specify_class.pointSource,
    #                      ['data_class', 'psf_class', 'lightModel', 'pointSource']]
    #     material = [fitting_specify_class.kwargs_data_joint['multi_band_list'],
    #                 fitting_specify_class.kwargs_model,
    #                 self.kwargs_result,
    #                 fitting_specify_class.kwargs_likelihood['image_likelihood_mask_list'],
    #                 kwargs_fixed_source, kwargs_fixed_ps, fitting_specify_class.kwargs_constraints,
    #                 fitting_specify_class.kwargs_numerics,
    #                 model_classes]
    #     pickle.dump([best_fit, chain_list_result, trans_paras, material], open(picklename, 'wb'))        
    
    def dump_result(self, tag = 'fitting_class'):
        dump_class = copy.deepcopy(self)
        if hasattr(dump_class.fitting_specify_class, 'data_process_class'):
            del dump_class.fitting_specify_class.data_process_class
        pickle.dump(dump_class, open(tag+'.pkl', 'wb'))        

        
def fitting_setting_temp(algorithm, fill_value_list = None):
    if algorithm == 'PSO':
        if fill_value_list is None:
            setting = {'sigma_scale': 0.8, 'n_particles': 150, 'n_iterations': 150}
        else:
            setting = {'sigma_scale': fill_value_list[0], 'n_particles': fill_value_list[1], 'n_iterations': fill_value_list[2]}
    elif algorithm == 'MCMC':     
        if fill_value_list is None:        
            setting = {'n_burn': 100, 'n_run': 200, 'walkerRatio': 10, 'sigma_scale': .1}
        else:
            setting = {'n_burn': fill_value_list[0], 'n_run': fill_value_list[1],
                       'walkerRatio': fill_value_list[2], 'sigma_scale': fill_value_list[3]}
    return setting