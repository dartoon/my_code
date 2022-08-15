#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 18:47:41 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from astropy.wcs import WCS
import glob
# from astropy.utils.data import get_pkg_data_filename

# folder = '/Users/Dartoon/Downloads/CEERS_JWST_data'
folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_JWST_MAST_data/'
all_files= glob.glob(folder+'/*clear*/*_i2d.fits')  #For NIRCam
# all_files= glob.glob(folder+'/jw*/*_i2d.fits')
filters = []
for file in all_files:
    try:
        filt = file.split('clear-')[1][:5]
    except:
        filt = file.split('_miri_')[1][:5]
    if filt not in filters:
        filters.append(filt)
    
#%% Make the plot
filt = 'f150w'
files = '/jw*/*{0}*_i2d.fits'.format(filt)
data_file = glob.glob(folder+files)
data_file.sort()
plt.figure(figsize=(12, 12))
fitsFile = pyfits.open(data_file[0])
fov_image = fitsFile[1].data # check the back grounp
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
wcs = WCS(header)
ax = plt.subplot(projection=wcs)
shape = np.shape(fov_image)
xs = [0, shape[1], shape[1],0 , 0]
ys = [0, 0, shape[0], shape[0], 0]
ax.plot(xs, ys, color="red", label = filt)
plt.text(np.mean(xs), np.mean(ys)+500, data_file[0].split('/')[-1].split('_')[1])
# =============================================================================
# # load the RA DEC for the second fits.
# =============================================================================
min_x, min_y = 0, 0
# for file in data_file:
for file in data_file[3:]:
    print(file)
    fitsFile = pyfits.open(file)
    fov_image = fitsFile[1].data # check the back grounp
    header_i = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    wcs_i = WCS(header_i)
    shape = np.shape(fov_image)
    positions_ = [(0,0), [shape[1],0], [shape[1],shape[0]], [0,shape[0]] ]
    sky = [wcs_i.pixel_to_world(positions_[i][0],positions_[i][1]) for i in range(4)]
    positions = np.array([wcs.world_to_pixel(sky[i]) for i in range(4)])
    xs_ = [positions[0][0], positions[1][0], positions[2][0], positions[3][0] , positions[0][0] ]
    ys_ = [positions[0][1], positions[1][1], positions[2][1], positions[3][1] , positions[0][1] ]
    min_x_ = positions[:,0].min()-100
    min_y_ = positions[:,1].min()-100
    if min_x > min_x_:
        min_x  = min_x_ 
    if min_y > min_y_:
        min_y  = min_y_ 
    ax.plot(xs_, ys_, color="red")
    plt.text(np.mean(xs_), np.mean(ys_), file.split('/')[-1].split('_')[1])
# =============================================================================
# # Plot the FOV for Second filter
# =============================================================================
filt = 'f356w'
files = '/jw*/*{0}*_i2d.fits'.format(filt)
data_file = glob.glob(folder+files)
data_file.sort()
for ct, file in enumerate(data_file):
    print(file)
    fitsFile = pyfits.open(file)
    fov_image = fitsFile[1].data # check the back grounp
    header_i = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    wcs_i = WCS(header_i)
    shape = np.shape(fov_image)
    positions_ = [(0,0), [shape[1],0], [shape[1],shape[0]], [0,shape[0]] ]
    sky = [wcs_i.pixel_to_world(positions_[i][0],positions_[i][1]) for i in range(4)]
    positions = np.array([wcs.world_to_pixel(sky[i]) for i in range(4)])
    xs_ = [positions[0][0], positions[1][0], positions[2][0], positions[3][0] , positions[0][0] ]
    ys_ = [positions[0][1], positions[1][1], positions[2][1], positions[3][1] , positions[0][1] ]
    min_x_ = positions[:,0].min()-100
    min_y_ = positions[:,1].min()-100
    label = None
    if ct == 0:
        label = filt
    ax.plot(xs_, ys_, color="green", label=label)
    plt.text(np.mean(xs_), np.mean(ys_), file.split('/')[-1].split('_')[1])

target_ID_list = []
target_ID_RA_DEC = []
target_ID_z = []

# =============================================================================
# #Plot target in fov
# =============================================================================
cata_file = '../../catalog_regions/SDSS_DR16Q_v4.fits'
hdul = pyfits.open(cata_file)
table = hdul[1].data
name = hdul[1].columns
table = table[table['RA']>213]
table = table[table['RA']<216]
table = table[table['Dec']>52]
table = table[table['Dec']<54]
targets = []
for i in range(len(all_files)):
    file = all_files[i]
    fitsFile = pyfits.open(file)
    fov_image = fitsFile[1].data # check the back grounp
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    wcs_i = WCS(header)
    for j in range(len(table)):
        RA, Dec, z = table[j][1], table[j][2],  table[j]['Z']
        pos = wcs_i.all_world2pix([[RA, Dec]], 1)[0]
        try:
            if pos[0]>0 and pos[1]>0 and fov_image[int(pos[1])-1, int(pos[0])-1 ] > 0 :
                # if z >2:#  or z_photo>2 and z_spec ==-99.0:
                pos_ = wcs.all_world2pix( [ [RA, Dec] ], 1)[0]
                target_ID = table[j][0]
                plt.scatter(pos_[0], pos_[1], c = 'red')
                plt.text(pos_[0], pos_[1], target_ID)
                if target_ID not in target_ID_list:
                    target_ID_list.append(target_ID)
                    target_ID_RA_DEC.append( [RA, Dec] )
                    target_ID_z.append([z, -99])
        except:
            None
f = open('../../catalog_regions/AEGIS_data_140612.csv',"r") ##This RA DEC of the optical counterparts is used to get hte reg file
string = f.read()
AEGIS_2014 = string.split('\n')   # Split in to \n
f = open('../../catalog_regions/AEGIS-XD_redshift_catalog.txt',"r") ##This RA DEC of the optical counterparts is used to get hte reg file
string = f.read()
AEGIS_redshift = string.split('\n')   # Split in to \n
def return_z(target_id):
    line = [AEGIS_redshift[i] for i in range(len(AEGIS_redshift)) if target_id in AEGIS_redshift[i]]
    if len(line)>1:
        print("find two same ID in redshift.")#, line)
    s = line[0]
    info = s.split(' ')
    info = [info_ for info_ in info if  info_ !='']
    if info[4] != '0' and float(info[3])>=0:
        z_spec = float(info[2])
    else:
        z_spec = -99.
    z_photo = float(info[6])
    return z_spec, z_photo
targets = []
for i in range(len(all_files)):
    file = all_files[i]
    fitsFile = pyfits.open(file)
    fov_image = fitsFile[1].data # check the back grounp
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    # plt_fits(fov_image)
    wcs_i = WCS(header)
    # print('search for file', file)
    for line in AEGIS_2014[17:-1]:
        info = line.split(' ')
        info = [info_ for info_ in info if  info_ !=''] 
        target_id, RA, Dec = info[0], float(info[3]), float(info[4])
        if RA != -99:
            pos = wcs_i.all_world2pix([[RA, Dec]], 1)[0]
        try:
            if pos[0]>0 and pos[1]>0 and fov_image[int(pos[1])-1, int(pos[0])-1 ] > 0 :
                z_spec, z_photo = return_z(target_id)
                
                if target_id not in target_ID_list:
                    target_ID_list.append(target_id)
                    target_ID_RA_DEC.append( [RA, Dec] )
                    target_ID_z.append([z_spec, z_photo])
                if z_spec >2:# or z_photo>2 and z_spec ==-99.0:
                    pos_ = wcs.all_world2pix( [ [RA, Dec] ], 1)[0]
                    plt.scatter(pos_[0], pos_[1], c = 'blue')  
                    plt.text(pos_[0], pos_[1], target_id)
                if z_spec <2 and z_photo>2:
                    pos_ = wcs.all_world2pix( [ [RA, Dec] ], 1)[0]
                    plt.scatter(pos_[0], pos_[1], c = 'green')  
                    plt.text(pos_[0], pos_[1], target_id)
                    # print(target_id, RA, Dec, 'flux:',fov_image[int(pos[1])-1, int(pos[0])-1 ], 'redshift:', z_spec, z_photo)
        except: 
            None

size = np.max(positions) * 6 
# plt.xlim(min_x,min_x+size)
# plt.ylim(min_y,min_y+size)
plt.tick_params(labelsize=25)
plt.legend(loc=1,prop={'size':15})
ax.set_xlabel('',fontsize = 25)
ax.set_ylabel('',fontsize = 25)
ax.set_aspect('equal', 'box')
plt.show()

# #%%
# write_file = open('target_info.txt','w') 
# write_file.write("#ID, RA, Dec, spec_z, photo_z\n")
# for i in range(len(target_ID_list)):
#     write_file.write("{0} {1} {2} {3} {4}\n".format(target_ID_list[i], target_ID_RA_DEC[i][0], target_ID_RA_DEC[i][1],
#                                                      target_ID_z[i][0], target_ID_z[i][1],))
# write_file.close()
