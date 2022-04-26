#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 15:39:02 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

filefolder = 'COSMOSweb_pointing/'
# filename = '1727_cosmos_mosaic_nircam_LW_exptime_scale1.0.fits'
filename = '1727_cosmos_mosaic_miri_exptime_scale1.0.fits'
fitsFile = pyfits.open(filefolder+filename)
reg_img = fitsFile[0].data # check the back grounp
header = fitsFile[0].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']

from astropy.wcs import WCS
wcs = WCS(header)

cata_folder = 'Chandra_COSMOS_Catalog/'
cata_file = 'chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits'
hdul = pyfits.open(cata_folder+cata_file)
table = hdul[1].data
name = hdul[1].columns

# for i in range(len(name)):  #Find that 24 and 25 are i-band countpart position
    # if 'z_best' in name[i].name:
    #     print(i)
#     if 'id_k_cfht' in name[i].name:
#         print(i, name[i].name, name[i+1], name[i+2], table[0][i+1], table[0][i+2])

#%%
cata_list = []
for i in range(len(table)):
    RA, Dec = table[i][24], table[i][25]
    if RA != -99:
        pos = wcs.all_world2pix([[RA, Dec]], 1)[0]
        try:
            exp = reg_img[int(pos[0]), int(pos[1]) ]
        except:
            exp = 0
        cata_list.append([i, pos[0], pos[1], exp, table[i][38], table[i][39]]) #idx, posx, posy, exp, z_best, z_spec
cata_list = np.array(cata_list)

plt.figure(figsize=(18,18))
plt.imshow(reg_img, origin='lower')
plt.scatter(cata_list[:,1], cata_list[:,2], c='r', s=10)
plt.show()

exps = cata_list[:,3]
z_bests = cata_list[:,4]
z_spec = cata_list[:,5]
print('numer of AGNs within', len(exps[exps!=0]))

plt.hist(exps[exps!=0])
plt.xlabel('exp time')
plt.show()

#%%
plt.hist(z_bests[ (exps!=0) * (z_bests>0)])
plt.xlabel('z_best')
plt.show()

plt.hist(z_spec[ (exps!=0) * (z_spec>0)])
plt.xlabel('z_spec')
plt.show()












