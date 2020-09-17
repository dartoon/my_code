#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 13:40:57 2020

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import copy
import pickle

class ResultPresent(object):
    """
    A class to quickly present the fitting result.
    """    
    def __init__(self, fitting_run_class, zp = 27.0, savename = 'result'):
        self.fitting_run_class = fitting_run_class
        self.delatPixel = fitting_run_class.fitting_specify_class.deltaPix
        self.zp = zp
        self.savename = savename
    
    def qso_final_plot(self, if_annuli=False, show_plot = False, arrows=False, save_plot = False):
        from decomprofile.tools_data.plot_tools import total_compare
        data = self.fitting_run_class.fitting_specify_class.kwargs_data['image_data']
        noise = self.fitting_run_class.fitting_specify_class.kwargs_data['noise_map']
        ps_list = self.fitting_run_class.image_ps_list
        ps_image = np.zeros_like(ps_list[0])
        for i in range(len(ps_list)):
            ps_image = ps_image+ps_list[i]
        galaxy_list = self.fitting_run_class.image_host_list
        galaxy_image = np.zeros_like(galaxy_list[0])
        for i in range(len(galaxy_list)):
            galaxy_image = galaxy_image+galaxy_list[i]
        model = ps_image + galaxy_image
        data_removePSF = data - ps_image
        norm_residual = (data - model)/noise
        flux_list_2d = [data, model, data_removePSF, norm_residual]
        label_list_2d = ['data', 'model', 'data-Point Source', 'normalized residual']
        
        flux_list_1d = [data, model, ps_image, galaxy_image]
        label_list_1d = ['data', 'model', 'QSOs', 'galaxy(s)']
        
        fig = total_compare(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d, delatPixel = self.delatPixel,
                      zp=self.zp, if_annuli=if_annuli, arrows= arrows, show_plot = show_plot)
        if show_plot == True:
            plt.show()
        else:
            plt.close()
            
        if save_plot == True:
            savename = self.savename
            fig.savefig(savename+"_qso_final_plot.pdf")
        
    def translate_result(self):
        import lenstronomy.Util.param_util as param_util
        from decomprofile.tools_data.measure_tools import model_flux_cal
        
        self.final_galaxy_result = copy.deepcopy(self.fitting_run_class.source_result)
        flux_sersic_model = model_flux_cal(self.final_galaxy_result)
        for i in range(len(self.final_galaxy_result)):
            source = self.final_galaxy_result[i]
            source['phi_G'], source['q'] = param_util.ellipticity2phi_q(source['e1'], source['e2'])
            source['flux_sersic_model'] = flux_sersic_model[i]
            source['flux_within_frame'] = np.sum(self.fitting_run_class.image_host_list[i])
            source['magnitude'] = -2.5*np.log10(source['flux_within_frame']) + self.zp
            self.final_galaxy_result[i] = source
         
        self.final_ps_result = copy.deepcopy(self.fitting_run_class.ps_result)
        for i in range(len(self.final_ps_result)):
            ps = self.final_ps_result[i]
            ps['flux_within_frame'] = np.sum(self.fitting_run_class.image_ps_list[i])
            ps['magnitude'] = -2.5*np.log10(ps['flux_within_frame']) + self.zp
        
    def dump_result(self):
        savename = self.savename
        dump_class = copy.deepcopy(self)
        if hasattr(dump_class.fitting_run_class.fitting_specify_class, 'data_process_class'):
            del dump_class.fitting_run_class.fitting_specify_class.data_process_class
        pickle.dump(dump_class, open(savename+'.pkl', 'wb'))      
            
