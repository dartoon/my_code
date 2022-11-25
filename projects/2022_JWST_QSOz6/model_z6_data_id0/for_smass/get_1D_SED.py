#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 23:48:59 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

target_id = 'J2255'

spec1d = pyfits.open('J2255_spec.fits')  
name_spec = spec1d[1].columns
unit =  spec1d[1].header['TUNIT3']
table_spec = spec1d[1].data
wave = table_spec['wave_model']
f_16 = table_spec['f_model_16']
f_50 = table_spec['f_model_50']
f_84 = table_spec['f_model_84']
plt.figure(figsize=(10, 6))

plt.plot(wave/10000., f_50, 'black', alpha=0.7)

plt.xlim(0.25, 6)
plt.ylim(0., 0.41)
    
plt.tick_params(labelsize=20)
plt.xlabel(r"$\lambda$ (um)",fontsize=27)
plt.ylabel(r"f$_\lambda$ 10$^{\rm" + " -{0}}}$ erg/s/cm$^2$/$\AA$".format(unit.split('e-')[1][:2]),fontsize=27)
plt.title(target_id,fontsize=27)

# plt.savefig('../figures/{0}_SED_map.pdf'.format(target_id[:5]))
# #plt.yticks([])


plt.show()