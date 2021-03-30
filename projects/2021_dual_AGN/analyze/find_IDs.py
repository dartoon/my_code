#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 18:08:19 2021

@author: Dartoon
"""

# import numpy as np
# import astropy.io.fits as pyfits
# import matplotlib.pyplot as plt
# from astropy.coordinates import SkyCoord
# from astropy import units as u
import glob

all_IDs  = glob.glob('../_John_fitted/*HSC-I/')
all_IDs = [all_IDs[i].split('/')[2].split('_HSC')[0] for i in range(len(all_IDs))]
all_IDs.sort()
# # fitted_IDs = glob.glob('../proofBHBH/model_Iband_zover_1/*')
# # fitted_IDs = [fitted_IDs[i].split('/')[-1] for i in range(len(fitted_IDs))]

# for ID in all_IDs: #Check if some bad ID exists.
    # # cha = ID.split('.')[-1]
    # cha = ID.split('.')[1].split('+')[0].split('-')[0]
    # if len(cha) > 2:
    #     print(ID)

from tools import read_info
for ID in all_IDs:
    # if ID not in fitted_IDs:
    RA, Dec, z = read_info(ID)
    # if RA == -99:
    #     pos = SkyCoord('{0} {1}'.format(ID[:2]+':'+ID[2:4]+':'+ID[4:9], ID[9:12]+':'+ID[12:14]+':'+ID[14:]), unit=(u.hourangle, u.deg))
    #     RA, Dec = pos.ra.degree, pos.dec.degree
    print(ID, RA, Dec, z)
    write_file = open('Detected_candidates_info.txt','r+')  #Will mv as material/ID_RA_DEC_z.txt
    write_file.read()
    write_file.write(ID+" {0} {1} {2}\n".format(RA, Dec, z))
    write_file.close()
