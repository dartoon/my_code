#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 21:20:09 2017

@author: dxh
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import astropy.io.fits as pyfits
import copy

from lenstronomy.Util import constants as const
import lenstronomy.Util.param_util as param_util

#file name:
filt='f160w'

import pickle
kwargs_lens_list, kwargs_lens_light_list, kwargs_source_list, kwargs_ps = pickle.load(open('sim_lens_noqso_ID_221/sim_kwargs.pkl','rb'))
