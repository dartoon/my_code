#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 12:17:37 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt


from astropy.wcs import WCS

def RA_Dec_in_fit(all_files, RA,Dec):
    filename = None
    for i in range(len(all_files)):
        file = all_files[i]
        print(file)
        fitsFile = pyfits.open(file)
        fov_image = fitsFile[1].data # check the back grounp
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        wcs_i = WCS(header)
        pos = wcs_i.all_world2pix([[RA, Dec]], 1)[0]
        print('pos')
        print(pos)
        try:
            if pos[0]>0 and pos[1]>0 and fov_image[int(pos[1])-1, int(pos[0])-1 ] > 0 :
                filename = file
        except:
            None
    return filename

def target_in_fits(all_files):
    target_info = []
    cata_file = '../catalog_regions/SDSS_DR16Q_v4.fits'
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
                    if z >2:
                        target_info.append([table[j][0], RA, Dec, [z], file])
            except:
                None
    f = open('../catalog_regions/AEGIS_data_140612.csv',"r") ##This RA DEC of the optical counterparts is used to get hte reg file
    string = f.read()
    AEGIS_2014 = string.split('\n')   # Split in to \n
    f = open('../catalog_regions/AEGIS-XD_redshift_catalog.txt',"r") ##This RA DEC of the optical counterparts is used to get hte reg file
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
        wcs_i = WCS(header)
        for line in AEGIS_2014[17:-1]:
            info = line.split(' ')
            info = [info_ for info_ in info if  info_ !=''] 
            target_id, RA, Dec = info[0], float(info[3]), float(info[4])
            if RA != -99:
                pos = wcs_i.all_world2pix([[RA, Dec]], 1)[0]
            try:
                if pos[0]>0 and pos[1]>0 and fov_image[int(pos[1])-1, int(pos[0])-1 ] > 0 :
                    z_spec, z_photo = return_z(target_id)
                    if z_spec >0 and z_spec>2:# or z_photo>2 and z_spec ==-99.0:
                        target_info.append([ target_id, RA, Dec, [z_spec, z_photo], file])
                    if z_spec <0 and z_photo>2:
                        target_info.append([ target_id, RA, Dec, [z_spec, z_photo], file])
            except: 
                None
    sdss_info = []
    aegis_zspec_info = []
    aegis_zphoto_info = []
    for i in range(len(target_info)):
        target_id, RA, Dec, z, file = target_info[i]
        if len(z) == 1:
            sdss_info.append(target_info[i])
        if len(z) == 2 and z[0]>2:
            aegis_zspec_info.append(target_info[i])
        if len(z) == 2 and z[0]<0:
            aegis_zphoto_info.append(target_info[i])
        
    return sdss_info + aegis_zspec_info + aegis_zphoto_info

