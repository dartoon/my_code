#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 13:23:13 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from functions_for_result import esti_smass, load_prop, load_info

# ID, mags, z = 'idx0', 
# 1,2,0,51,35
idx = 1
# root_folder = '../*/*'  #Include HST
root_folder = './*'
mag_result = load_prop(idx, root_folder = root_folder)
target_id, z = load_info(idx)
#%%
# esti_smass(ID = '202208'+str(idx), mags_dict = mag_result, z = z, flag = 1, if_run_gsf=True)
print('Run estimate')
import time
t1 = time.time()
esti_smass(ID = '100'+str(idx), mags_dict = mag_result, z = z, flag = 1, if_run_gsf=True)
t2 = time.time()