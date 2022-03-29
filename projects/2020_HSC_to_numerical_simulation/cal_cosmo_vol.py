#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:25:18 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import numpy as np
# from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

zs = 0.5
dis = cosmo.comoving_distance(np.array([zs]))  

# HSC_cov = 1000 #sq deg
HSC_cov = 500 * 3600**2 #sq arcsec 

dia = cosmo.angular_diameter_distance(zs)
dia_sq = dia ** 2 

vol = dia_sq * dis[0]
print(vol.value)