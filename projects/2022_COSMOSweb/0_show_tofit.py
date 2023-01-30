#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 22:50:42 2023

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import warnings, pickle
warnings.filterwarnings("ignore")
from astropy.wcs import WCS
from galight.tools.astro_tools import plt_many_fits
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import pickle, copy, glob
import time
filt_i = 1
filt = ['F115W', 'F150W','F277W', 'F444W'][filt_i]
cata_list = pickle.load(open('material/cata_list.pkl','rb'))

not_fit_id = []
for i in range(0, len(cata_list)):
# for i in [0]:
    savename = 'fit_material/fit_notrunyet_{2}_idx{0}_psf{1}.pkl'.format(i,4,filt)
    try:
        fit_run = pickle.load(open(savename,'rb'))
        print('idx {0} for fitting:'.format(i))
        f = fit_run.fitting_specify_class.plot_fitting_sets()
    except:
        print('idx {0} is skiped'.format(i))
        not_fit_id.append(i)
    time.sleep(0.5)
    
    # # Quick run:
    # fit_run.run(algorithm_list = ['PSO'], fitting_level=['norm'])
    # fit_run.plot_final_qso_fit(target_ID =None)