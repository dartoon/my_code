#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 16:44:10 2022

@author: Dartoon
"""

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

reg_file = 'catalog_regions/AEGIS-XD_optical_pix.reg'
f = open(reg_file,"r")
string = f.read()
lines = string.split('\n')   # Split in to \n


positions = []
for i in range(len(lines)):
    if 'circle' in lines[i]:
        x = lines[i].split('circle(')[1].split(',')[0]
        y = lines[i].split('circle(')[1].split(',')[1]
        positions.append([float(x),float(y)])
positions = np.array(positions)

in_fov_idx = []
for i in range(len(positions)):
    pos = positions[i]
    try:
        exp = reg_img[int(pos[1]), int(pos[0]) ]
    except:
        exp = 0
    if exp!=0:
        in_fov_idx.append(i)  #in_fov_idx record the id that have exp >0

#In the following, we make cross match to get catalog information.
reg_file_deg = 'catalog_regions/AEGIS-XD_optical.reg' #This is the origial RA Dec from the AEGIS_data_140612.csv file.
f = open(reg_file_deg,"r")
string = f.read()
lines_deg = string.split('\n')   # Split in to \n
        
pos_deg = []
for i in range(len(lines_deg)):
    if 'circle' in lines_deg[i]:    
        ra = lines_deg[i].split('circle(')[1].split(',')[0]
        dec = lines_deg[i].split('circle(')[1].split(',')[1]
        pos_deg.append([ra, dec])


f = open('catalog_regions/AEGIS_data_140612.csv',"r") ##This RA DEC of the optical counterparts is used to get hte reg file
string = f.read()
AEGIS_2014 = string.split('\n')   # Split in to \n

in_aegids = []
for i in range(len(in_fov_idx)):
    for j in range(len(AEGIS_2014)):
        if pos_deg[in_fov_idx[i]][0] in AEGIS_2014[j] and pos_deg[in_fov_idx[i]][1] in AEGIS_2014[j]:
            in_aegids.append(AEGIS_2014[j].split(' ')[0])

f = open('catalog_regions/AEGIS-XD_redshift_catalog.txt',"r") ##This RA DEC of the optical counterparts is used to get hte reg file
string = f.read()
AEGIS_redshift = string.split('\n')   # Split in to \n
# #%% Check if the number of sample with spec-z is actually be around 365... Yes, it is 365
# c = 0
# for s in AEGIS_redshift:
#     if 'aegis_' in s:
#         info = s.split(' ')
#         info = [info_ for info_ in info if  info_ !='']
#         if info[4] != '0' and float(info[3])>=3:
#             print(info[0], info[2])
#             c = c+1
#print(c)
for _id in in_aegids:
    for s in AEGIS_redshift:
        if _id in s:
            info = s.split(' ')
            info = [info_ for info_ in info if  info_ !='']
            if info[4] != '0' and float(info[3])>=0 and float(info[2])>2:
                i = [i for i in range(len(in_aegids)) if _id == in_aegids[i]][0]
                # pos_deg[in_fov_idx[i]]
                print(_id, info[2], pos_deg[in_fov_idx[i]])
            # print("Let's check photo z")
            elif float(info[6])>2: #Check up for photo-z.
                i = [i for i in range(len(in_aegids)) if _id == in_aegids[i]][0]
                # pos_deg[in_fov_idx[i]]
                print(_id, info[6], info[2], info[3], pos_deg[in_fov_idx[i]],'pz')
                
# #%% Check how many photo-z (only) over 2:
# count = 0
# count1  = 0
# for s in AEGIS_redshift:
#     if 'aegis_' in s:
#         count = count +1
#         info = s.split(' ')
#         info = [info_ for info_ in info if  info_ !='']
#         if float(info[2])>2:
#             count1 = count1+1
#             # i = [i for i in range(len(in_aegids)) if info[0] == in_aegids[i]][0]
#             # pos_deg[in_fov_idx[i]]
#             print(info[0], info[2])               
                
                
    
