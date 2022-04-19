#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 17:31:22 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pandas as pd
filename = 'JWST_CEERS/ceers5_f150w_cat.ecsv'

sample = pd.read_csv(filename)

# ID_list = sample['object_id_1']
# RA_list = sample['ra_1']
# Dec_list = sample['dec_1']