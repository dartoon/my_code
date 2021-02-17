#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 14:27:40 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy import units as u


f = open('data/Reines_2015_table_1.txt',"r")
# ID  Running identification number
# NSAID NASA-Sloan Atlas identification number
# SDSS SDSS identification
# P-MJD-F Plate-MJD-Fiber of SDSS spectrum
# z  Redshift (1)
# mag  iMag Absolute i-band magnitude (2)
# mag  g-i  The (g-i) color (2)
# [solMass] logM* Log stellar mass; corrected for AGN contribution (3)
# [solMass] logMBH Log black hole mass (3)
string = f.read()
Reines_t1 = string.split('\n')   # Split in to \n
Reines_t1 = [Reines_t1[i].split(' ') for i in range(len(Reines_t1)) if Reines_t1[i][0]!= '#']

for i in range(len(Reines_t1)):
    ID = Reines_t1[i][2]
    Ra = ID[1:10]
    Ra = Ra[:2] + ' ' + Ra[2:4] + ' ' + Ra[4:]
    Dec = ID[10:]
    Dec = Dec[:3] +  ' ' + Dec[3:5] + ' ' + Dec[5:]
    pos = SkyCoord('{0}{1}'.format(Ra, Dec), unit=(u.hourangle, u.deg))
    RA, DEC = pos.ra.degree, pos.dec.degree
    print(ID, round(RA,6), round(DEC,6) )