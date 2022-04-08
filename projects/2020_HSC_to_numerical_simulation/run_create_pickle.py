#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 16:22:19 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from sum_prep import load_HSC_comp, load_HST_comp, imf_dict
import pickle

# Imag_list = {0.3: 20.0, 0.5: 20.5, 0.7: 21.5}

# zs = 1.5
# obs_dict, sim_dict = load_HST_comp(no_noise=False)
# pickle.dump([obs_dict, sim_dict], open('pickles/comb_zs{0}.pkl'.format(zs), 'wb'))  
# obs_dict, sim_dict = load_HST_comp(no_noise=True)
# pickle.dump([obs_dict, sim_dict], open('pickles/comb_zs{0}_nonoise.pkl'.format(zs), 'wb'))  


## For ten time realiszation.
# for i in range(10):
#     for zs in [0.3, 0.5, 0.7]:
#         if zs<1:
#             Imag= Imag_list[zs]
#             obs_dict, sim_dict = load_HSC_comp(zs =zs, I_mag_break=Imag, no_noise=False)
#             pickle.dump([obs_dict, sim_dict], open('pickles/comb_zs{0}_Imag_{1}_{2}.pkl'.format(zs,Imag,i), 'wb'))    

# #No noise
# for zs in [0.7]:
#     Imag= Imag_list[zs]
#     obs_dict, sim_dict = load_HSC_comp(zs =zs, I_mag_break=Imag, no_noise=True)
#     pickle.dump([obs_dict, sim_dict], open('pickles/comb_zs{0}_Imag_{1}_nonoise.pkl'.format(zs,Imag), 'wb'))    

# #For Imag = 21.0 for MBII at 0.6
# zs = 0.5
# Imag = 21.0
# obs_dict, sim_dict = load_HSC_comp(zs =zs, I_mag_break=Imag, no_noise=True)
# pickle.dump([obs_dict, sim_dict], open('pickles/comb_zs{0}_Imag_{1}_nonoise.pkl'.format(zs,Imag), 'wb'))    
# obs_dict, sim_dict = load_HSC_comp(zs =zs, I_mag_break=Imag, no_noise=False)
# pickle.dump([obs_dict, sim_dict], open('pickles/comb_zs{0}_Imag_{1}.pkl'.format(zs,Imag), 'wb'))    


# #For Imag = 21.5 for Imag different comparison 
zs = 0.5
Imag = 22.0
obs_dict, sim_dict = load_HSC_comp(zs =zs, I_mag_break=Imag, no_noise=True)
pickle.dump([obs_dict, sim_dict], open('pickles/comb_zs{0}_Imag_{1}_nonoise.pkl'.format(zs,Imag), 'wb'))    
obs_dict, sim_dict = load_HSC_comp(zs =zs, I_mag_break=Imag, no_noise=False)
pickle.dump([obs_dict, sim_dict], open('pickles/comb_zs{0}_Imag_{1}.pkl'.format(zs,Imag), 'wb'))    
