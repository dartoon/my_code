#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 20:51:21 2018

@author: dartoon

Comparing the difference between difference stars.
"""

import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pylab as plt

O_star = pyfits.getdata('O_star_test00.fits')
K_star = pyfits.getdata('K_star_test00.fits')

#plt.imshow(O_star, origin='low')
plt.imshow((K_star-O_star)/K_star, origin='low', vmax=0.4, vmin=-0.4)
plt.colorbar()
plt.show()

print "The precent of mismatch:", round(np.sum(abs(K_star-O_star))/np.sum(K_star),3)