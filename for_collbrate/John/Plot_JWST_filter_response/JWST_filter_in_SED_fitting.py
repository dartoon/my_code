#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 00:03:18 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'


plt.figure(figsize=(10, 6))
filt_id = {'F115W': '352', 'F150W': '353', 'F200W': '354', 
           'F277W': '355', 'F356W': '356', 'F444W': '357', 'F410M': '362'}
ivd = {v: k for k, v in filt_id.items()}
plt.xlim(0.25, 6)
plt.ylim(0., 10)
xmin, xmax, ymin, ymax = plt.axis()

for i, fid in enumerate(ivd.keys()):
    f_fil = np.loadtxt('./{0}.fil'.format(fid))
    f_fil[:,2] = f_fil[:,2]/f_fil[:,2].max() * (ymax-ymin) * 0.1
    plt.plot(f_fil[1:,1]/10000., f_fil[1:,2], label='{0}'.format(ivd[fid]))
plt.legend(prop={'size':15}, ncol=2, loc = 1)
plt.tick_params(labelsize=20)
plt.xlabel(r"$\lambda$ (um)",fontsize=27)
plt.ylabel(r"f$_\lambda$ 10$^{\rm" + " -{0}}}$ erg/s/cm$^2$/$\AA$".format(19),fontsize=27)
plt.show()


