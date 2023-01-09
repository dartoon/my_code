#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 09:39:13 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from gsf import gsf

folder = 'esti_smass/202212311012' 
gsf.run_gsf_all(folder+'/sample.input', 6, idman=None)