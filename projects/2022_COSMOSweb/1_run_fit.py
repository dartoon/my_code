#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 11:32:37 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle,glob

# import sys
# count = int(sys.argv[1]) - 1 # 1 - 12586
count = 1

# filt_i = 1
# filt = ['F115W', 'F150W','F277W', 'F444W'][filt_i]
cata_list = pickle.load(open('material/cata_list.pkl','rb'))

# fit_files = glob.glob('fit_material/fit_notrunyet_*.pkl')
fit_files = glob.glob('fit_material/fit_notrunyet_F115W_idx10*.pkl')
fit_files.sort()

# i = count
for i in range(len(fit_files)):
    fit_run = pickle.load(open(fit_files[i],'rb'))
    # Quick run:
    fit_run.run(algorithm_list = ['PSO','PSO','PSO'], fitting_level=['norm', 'deep','deep'])
    fit_run.plot_final_qso_fit(target_ID =None, show_plot=False)
    savename = fit_files[i].replace('_notrunyet_', '_run_')
    savename = savename.replace('fit_material', 'fit_result')
    pickle.dump(fit_run , open(savename, 'wb'))