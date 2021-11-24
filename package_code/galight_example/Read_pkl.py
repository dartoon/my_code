#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 19:49:11 2020

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits


# Test load pkl
import pickle
picklename = 'savename.pkl'
fitting_run_class = pickle.load(open(picklename,'rb'))
fitting_run_class.plot_all()
fitting_run_class.run_diag()

