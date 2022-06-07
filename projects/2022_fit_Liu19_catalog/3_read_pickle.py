#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 14:22:19 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
from galight.tools.astro_tools import plt_fits
ID = 10004
band = 'I'
folder = 'fit_result/'
file_ = glob.glob(folder+"{0}-{1}.pkl".format(ID, band))

if file_ != []:
    file = file_[0]
    fit_run = pickle.load(open(file,'rb'))
    print(fit_run.final_result_galaxy)
    host_image = fit_run.flux_2d_out['data'] - fit_run.image_ps_list[0]
    plt_fits(host_image)
    