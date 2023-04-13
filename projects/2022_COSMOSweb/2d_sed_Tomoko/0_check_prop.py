#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 10:28:04 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

f = open('fmos_alma_cosmosweb.cat','r')
string = f.read()
lines = string.split('\n')
lines = [lines[i] for i in range(len(lines)) if 'FMOS_J09' in lines[i]]

