#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 16:05:10 2020

@author: dartoon
"""

import numpy as np
line_means = ['id', 'z', 'ra', 'dec', 'fix_sersic_n', 'sersic_n_fitted', 'sersic_re_fitted', 'sersic_n_corrected',
         'sersic_re_corrected', 'host_mag_g', 'host_mag_r', 'host_mag_i', 'host_mag_z', 'host_mag_y',
         'ps_mag_g', 'ps_mag_r', 'ps_mag_i', 'ps_mag_z', 'ps_mag_y', 'decomposition_chisq', 'stellar_mass', 
         'sed_chisq', 'logMBH', 'logMBH_err']
infers  = np.loadtxt('./sdss_quasar_decomposition_v1.txt', dtype=str)
HSC_z = infers[:,1].astype(np.float)
HSC_Mstar = infers[:,20].astype(np.float)
HSC_MBHs = infers[:,22].astype(np.float)
HSC_MBHs_err = infers[:,23].astype(np.float)

