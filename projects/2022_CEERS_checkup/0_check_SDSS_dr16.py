#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:22:30 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt


filefolder = 'JuneMaps/'
filename = 'CEERS_nircam_lw_June2022.fits'
fitsFile = pyfits.open(filefolder+filename)
reg_img = fitsFile[0].data # check the back grounp
header = fitsFile[0].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
# from astropy.wcs import WCS
# wcs = WCS(header)

cata_file = 'catalog_regions/SDSS_DR16Q_v4.fits'
hdul = pyfits.open(cata_file)
table = hdul[1].data
name = hdul[1].columns

table = table[table['RA']>213]
table = table[table['RA']<216]
table = table[table['Dec']>52]
table = table[table['Dec']<54]

# write_file = open('SDSS_dr16.reg','w') 
# for i in range(len(table)):
#     write_file.write("circle({0},{1},0.005)\n".format(table[i]['RA'],table[i]['Dec']) )
# write_file.close()

# from astropy.coordinates import SkyCoord
# import astropy.units as u
# from regions import CircleSkyRegion
# from astropy.coordinates import Angle, SkyCoord
# for i in range(len(table)):
#     RA, Dec = table[i][1], table[i][2]
#     center = SkyCoord(RA,Dec, unit='deg')
#     radius = Angle(0.001, 'deg')
#     region = CircleSkyRegion(center, radius)
#     region.write('SDSS_regions/my_region{0}.reg'.format(i), overwrite=True)


# Find out the time.
# # ra, dec = 215.0233234, 53.01020815  #z = 1.646
# # ra, dec = 215.0359019, 53.00112292 # z = 2.588
# ra, dec = 214.9316085, 52.90871718  # z = 3.442
# table = hdul[1].data
# table = table[table['RA']>ra-0.001]
# table = table[table['RA']<ra+0.001]
# table = table[table['Dec']>dec-0.001]
# table = table[table['Dec']<dec+0.001]
# print(table[0]['Z'])

#%%
# import scipy.optimize as op
# # RA, Dec = (215.0233234, 53.01020815)
# RA, Dec = (215.0359019, 53.00112292 )
# # RA, Dec = (214.9316085, 52.90871718  )
# def cal_dis(theta,RA,Dec,wcs):
#     x, y = theta
#     cal_RA,cal_Dec = wcs.pixel_to_world_values(x,y)
#     dis = np.sqrt( (Dec - cal_Dec)**2 + (RA-cal_RA)**2 ) 
#     return dis

# def find_xy(RA, Dec, wcs, x=23310.317, y=7331.0882):
#     """
#     Use minimazation to find the position to match RA Dec. 
#     """
#     theta = x, y
#     result = op.minimize(cal_dis, theta, args=(RA,Dec,wcs), method='nelder-mead',
#             options={'xatol': 1e-1, 'disp': False})
#     x, y = result['x']
#     return x, y
# x, y = find_xy(RA, Dec, wcs)
# print(x,y)
# print(wcs.pixel_to_world_values(x,y), RA, Dec)
# print(cal_dis([x,y], RA, Dec, wcs))
# print(reg_img[int(y), int(x)])

# #%%
# from astropy.coordinates import SkyCoord
# import astropy.units as u
# from regions import CircleSkyRegion
# from astropy.coordinates import Angle, SkyCoord


# cata_list = []
# for i in range(len(table)):
#     RA, Dec = table[i][1], table[i][2]
#     if RA != -99:
#         # pos = wcs.all_world2pix([[RA, Dec]], 1)[0]
#         pos = find_xy(RA,Dec,wcs)
#         if pos[0]>0 and pos[1]>0:
#             try:
#                 exp = reg_img[int(pos[1]), int(pos[0]) ]
#             except:
#                 exp = 0
#         else:
#                 exp = 0
#         if exp!=0:
#             center = SkyCoord(RA,Dec, unit='deg')
#             radius = Angle(0.001, 'deg')
#             region = CircleSkyRegion(center, radius)
#             region.write('SDSS_regions/my_region{0}.reg'.format(i), overwrite=True)
            
#         cata_list.append([i, pos[0], pos[1], exp, table[i]['Z']]) #idx, posx, posy, exp, z_best, z_spec
# cata_list = np.array(cata_list)
# exps = cata_list[:,3]
# z_bests = cata_list[:,-1]

# plt.figure(figsize=(18,18))
# plt.imshow(reg_img, origin='lower')
# plt.scatter(cata_list[:,1][exps!=0], cata_list[:,2][exps!=0], c='r', s=10)
# plt.show()

# print('numer of AGNs within', len(exps[exps!=0]))
# plt.hist(exps[exps!=0])
# plt.xlabel('exp time')
# plt.show()

# #%%
# plt.hist(z_bests[ (exps!=0) * (z_bests>0)])
# plt.xlabel('z_best')
# plt.show()


