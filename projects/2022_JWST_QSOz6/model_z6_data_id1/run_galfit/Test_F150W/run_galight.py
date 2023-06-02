#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:11:34 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from galight.tools.astro_tools import read_pixel_scale
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess

filt = 'F150W'


for i in range(1,6):
    fitsFile = pyfits.open('J2236_{0}_data_cutout.fits'.format(filt))
    #Load the QSO image data:
    qso_image = fitsFile[0].data #target is at the image center of the cutout
    #Load header for the zp information
    header = fitsFile[0].header 
    #Load noise level map:
    fitsFile = pyfits.open('J2236_{0}_noise.fits'.format(filt))
    data_error= fitsFile[0].data 
    #Calculate the zeropoint
    pixscale = read_pixel_scale(header)  #Derive pixel scale
    zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631)  #For AB magnitude
    print('pixel scale', pixscale, 'ABzeropoint', zp)
    
    #Load the PSF data: 
    PSF = pyfits.getdata('PSF_{0}_top{1}.fits'.format(filt,i)) 
    # In[4]:
    #Prepare the fitting material
    data_process = DataProcess(fov_image = qso_image, fov_noise_map = data_error, 
                               target_pos = [len(qso_image)/2]*2, header = header,
                               rm_bkglight = False, if_plot=False, zp = zp)
    data_process.generate_target_materials(radius='nocut')
    
    #Manually input the PSF:
    data_process.PSF_list = [PSF]
    
    print("Show the flux profile (1D) of our QSO and the adopted PSF.")
    # Compare the 1D profile of all the components.
    data_process.profiles_compare(norm_pix = 5, if_annuli=False, y_log = False,
                      prf_name_list = (['target'] + ['PSF{0}'.format(i) for i in range(len(data_process.PSF_list))]) )
    
    #Check if all the materials is given, if so to pass to the next step.
    data_process.checkout() #Check if all the materials is known.
    
    
    # ## In the following settings, we fix the Sersic_n = 1 to perform the fitting.
    
    # In[5]:
    
    
    #Pass the data_process to FittingSpeficy
    fit_sepc = FittingSpecify(data_process)
    #Prepare the fitting sequence, keywords see notes above.
    fit_sepc.prepare_fitting_seq(point_source_num = 1, fix_n_list= [[0,1]], supersampling_factor = 5)  #i.e., the comp0's Sersic n is fixed to 1.
    ##One can also test to use the following line to set a free sersic n to fit
    #fit_sepc.prepare_fitting_seq(point_source_num = 1) 
    fit_sepc.kwargs_numerics['point_source_supersampling_factor'] = 1
    
    #Plot the initial settings for fittings. 
    fit_sepc.plot_fitting_sets()
    
    #Build up and to pass to the next step.
    fit_sepc.build_fitting_seq()
    
    
    # In[6]:
    
    fit_run = FittingProcess(fit_sepc, savename = '{0}_fit_psf{1}_ssf1'.format(filt,i), fitting_level=['norm','deep']) 
    fit_run.run(algorithm_list = ['PSO', 'PSO'], setting_list = None)
    fit_run.plot_final_qso_fit()
    fit_run.dump_result()
    
    # In[8]:
    print(filt, fit_run.final_result_galaxy[0])
    print(filt, fit_run.final_result_ps[0]['magnitude'])
    print(fit_run.reduced_Chisq)
    print(fit_run.zp)
    
    
    # In[8]:
    
    
    print( -2.5*np.log10(np.sum(fit_run.flux_2d_out['data'] - fit_run.image_ps_list[0]) ) + zp)
