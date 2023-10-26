#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 21:48:42 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from galight.tools.plot_tools import profile_plots

import pickle
import glob

# name = 'dual_result-band-I-s21a'
# picklename = 'fit_result/{0}.pkl'.format(name)
picklename = glob.glob('fit_result/*.pkl')[3]
name = picklename.split('/')[1].split('.pkl')[0]
fit_run = pickle.load(open(picklename,'rb'))
fit_run.plot_final_galaxy_fit()
fit_run.fitting_specify_class.plot_fitting_sets()
host = fit_run.image_host_list
AGN = fit_run.image_ps_list
# bulge_Re1 = fit_run.final_result_galaxy[0]['R_sersic']
# disk_Re1 = fit_run.final_result_galaxy[1]['R_sersic']
flux_list_2d = host + AGN
label_list_2d = ['host{0}'.format(i) for i in range(len(host))] + ['AGN{0}'.format(i) for i in range(len(AGN))]
flux_list_1d = flux_list_2d 
label_list_1d = label_list_2d
profile_plots(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d,
              deltaPix = fit_run.fitting_specify_class.deltaPix,
              target_ID =  picklename, if_annuli=True)   