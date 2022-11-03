#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 17:23:49 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt


import glob, pickle
run_folder = 'stage3_all*/' #!!!
filt = 'F356W'
idx = 0
PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
from galight.tools.astro_tools import plt_fits

# PSF_lib_files = glob.glob('../../2022_CEERS_checkup/real_analysis/model_JWST_stage3_Masafusa/material/*'+filt[:-1]+'*_PSF_Library.pkl')[0]
# PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))


for i in range(len(PSF_list_clean)):
    plt_fits(PSF_list_clean[i][40:-40, 40:-40])