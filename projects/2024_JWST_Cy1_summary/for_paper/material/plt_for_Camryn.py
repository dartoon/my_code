#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 17:11:03 2024

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy, matplotlib
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')
from astropy.coordinates import Angle
from photutils.aperture import RectangularAperture
from matplotlib.colors import LogNorm

target_info = {'0': {'target_id': 'J2255+0251', 
                      'RA': '22:55:38.04',
                      'Dec': '+02:51:26.6',
                      'z': 6.34,
                      'theta': 137.475, #The NIRSpec filter angle
                      'M1450': -23.9
                      } ,
               '1': {'target_id': 'J2236+0032', 
                        'RA': 339.1857947085765, #'22:36:44.5800',
                        'Dec':  0.5490648434896199, #'+00 32 56.90',
                        'z': 6.40,
                        'theta': 139.818,
                        'M1450': -23.8
                        } ,
               '2': {'target_id': 'J0844-0132',   #It's actually J0844-0132 but first use as J0844-1321
                     'RA': '08:44:08.61',
                     'Dec':  '-01:32:16.5',
                     'z': 6.18,
                     'theta': 138.521,
                      'M1450': -23.7
                     } ,
               '3': {'target_id': 'J0911+0152', 
                     'RA': '09:11:14.27',
                     'Dec':  '+01:52:19.4',
                     'z': 6.07,
                     'theta': None,
                      'M1450': -22.1
                     } ,
               '4': {'target_id': 'J0918+0139', 
                     'RA': '09:18:33.17',
                     'Dec':  '+01:39:23.4',
                     'z': 6.19,
                     'theta': 138.542,
                      'M1450': -23.7
                     } ,
               '5': {'target_id': 'J1425-0015', 
                     'RA': '14:25:17.705',
                     'Dec':  '-00:15:40.87',
                     'z': 6.18,
                     'theta': 148.04,
                      'M1450': -23.4} ,
               '6': {'target_id': 'J1512+4422', 
                     'RA': '15:12:48.7100',
                     'Dec':  '+44:22:17.5',
                     'z': 6.18,
                     'theta': 171.636,
                      'M1450': -23.1
                     } ,
               '7': {'target_id': 'J1525+4303', 
                     'RA': '15:25:55.7900 ',
                     'Dec':  '+43:03:24.0',
                     'z': 6.27,
                     'theta': 138.76,
                      'M1450': -23.9} ,
               '8': {'target_id': 'J1146-0005', 
                     'RA': '11:46:58.89',
                     'Dec':  '-00:05:37.7',
                     'z': 6.3,
                     'theta': None,
                      'M1450': -21.5
                     } ,
               '9': {'target_id': 'J1146+0124', 
                     'RA': '11:46:48.42',
                     'Dec':  '+01:24:20.1',
                     'z': 6.27,
                     'theta': None,
                      'M1450': -23.7
                     } ,
                      } 

idx = 9  #Change idx from 0 to 9 
fig, ax = plt.subplots( figsize=(5, 5))
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
image_file = 'data_minus_QSO/{0}_host.fits'.format(target_id)
image = pyfits.getdata(image_file)
deltaPix = 0.03149
im1 = ax.imshow(image, origin='lower', norm=LogNorm(vmin=0.0001,vmax=3),cmap = my_cmap) 
if target_info[str(idx)]['theta'] != None:
    theta = Angle(target_info[str(idx)]['theta'], 'deg')
    f_center = len(image)/2
    w = 0.2 / deltaPix
    h = 0.6 / deltaPix
    aper = RectangularAperture((len(image)/2, len(image)/2), w, h, theta=theta)
    aper.plot(color='white',
              lw=0.8,axes=ax)
ax.set_title(target_id, fontsize=17, fontweight="bold")  # Add title
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False) 
