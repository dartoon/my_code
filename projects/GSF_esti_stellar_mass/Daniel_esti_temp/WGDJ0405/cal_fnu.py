#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 23:12:18 2020

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits


mag = [20.2, 22.0, 25.0]
mag_err = [0.1, 0.1, 0.3]

fnu = [10 ** ((mag[i]-25)/(-2.5)) for i in range(len(mag))]
fnu_up = [10 ** ((mag[i]-mag_err[i]-25)/(-2.5)) for i in range(len(mag))]
fnu_dw = [10 ** ((mag[i]+mag_err[i]-25)/(-2.5)) for i in range(len(mag))]
fnu_err = [(fnu_up[i]-fnu_dw[i])/2 for i in range(len(mag))]

print([(round(fnu[i],3), round(fnu_err[i],3)) for i in range(len(mag))])