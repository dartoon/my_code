#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 20:30:48 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import os
import glob

# fit_run = pickle.load(open('fit_material/fit_run_idx9_.pkl','rb'))
# fit_run.plot_final_qso_fit()
# f = open("../model_JWST_stage3_Masafusa/target_idx_info.txt","r")
# string = f.read()
# lines = string.split('\n')   # Split in to \n

dp_files = glob.glob('fit_material/data_process_idx*.pkl')
dp_files.sort()
for i in range(len(dp_files)):
# for i in range(0,100):
    file = glob.glob('fit_material/'+'fit_run_idx*_*_psf*_{0}.pkl'.format(i))[0]
    if file == []:
        print(i)
    # # print(dp_files[i], file)
    # new_filename = dp_files[i].replace('data_process', 'fit_run')[:-4]+'_{0}.pkl'.format(i)
    # # print(file, new_filename)
    # os.rename(file, new_filename)


    