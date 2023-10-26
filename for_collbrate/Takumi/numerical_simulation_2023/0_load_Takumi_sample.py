#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 13:55:12 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('matched.csv')

df.keys()

MBH = df['logMBH_1']
Mass = df['M_st5']
logLbol = df['logLbol_final']

    
logLedd_overall = 38. + np.log10(1.2) + MBH
Eddr_overall = logLbol - logLedd_overall

plt.scatter(Mass, MBH, c = logLbol)
plt.xlim([9,13])
plt.ylim([6,10.5])
plt.show()

plt.scatter(MBH, Eddr_overall, c = logLbol)
plt.show()
