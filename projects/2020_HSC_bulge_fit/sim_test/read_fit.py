#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 14:49:45 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle

import glob
files = glob.glob("sim_result_bulge_n4/round0_ID9_*.pkl")
result = pickle.load(open(files[0],'rb'))

