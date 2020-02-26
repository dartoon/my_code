#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 22:49:33 2020

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import os
#import numpy as np
from mirage.catalogs import catalog_generator
from mirage.catalogs import create_catalog

ra = np.random.random(5) + 80.
dec = np.random.random(5) - 69.8
mags1 = np.random.random(5) + 17.
mags2 = np.random.random(5) + 18.5

ptsrc = catalog_generator.PointSourceCatalog(ra=ra, dec=dec)
print('RA: {}'.format(ptsrc.ra))
print('Dec: {}'.format(ptsrc.dec))

ptsrc.add_magnitude_column(mags1, instrument='nircam', filter_name='f090w', magnitude_system='abmag')
ptsrc.add_magnitude_column(mags1+0.1, instrument='nircam', filter_name='f444w', magnitude_system='abmag')
ptsrc.add_magnitude_column(mags2, instrument='niriss', filter_name='f200w', magnitude_system='abmag')

ptsrc.table
ptsrc.save('ptsrc_test.cat')
ptsrc.location_units
x_pix = np.random.random(5) * 2048
y_pix = np.random.random(5) * 2048
ptsrc_pix = catalog_generator.PointSourceCatalog(x=x_pix, y=y_pix)
ptsrc_pix.add_magnitude_column(mags1, instrument='nircam', filter_name='f090w')

ptsrc_pix.location_units
ptsrc_pix.table

#Galaxy:

radius = np.random.random(5) + 0.5
ellip = np.random.random(5) + 0.45
posang = np.random.random(5) + 27.
sersic = np.random.random(5) + 3.3

gal = catalog_generator.GalaxyCatalog(ra=ra, dec=dec, ellipticity=ellip, radius=radius, sersic_index=sersic,
                                      position_angle=posang, radius_units='arcsec')

gal.add_magnitude_column(mags1, instrument='nircam', filter_name='f444w', magnitude_system='stmag')
gal.table
gal.save('galaxy_test.cat')