#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 17:58:21 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from ID_list import ID_list
import glob

from ast import literal_eval
def read_string_list(string):
    """
    Translate a string-list-dict to a list-dict. 
    Not allow array inside the dict...
    """
    string = ''.join(string.split('array(['))
    string = ''.join(string.split('])'))    
    string = string.replace('nan', '-99')
    string = string.replace('inf', '-99')
    string_list = string.split('{')[1:]
    string_list = [literal_eval("{"+string_list[i].split('}')[0]+"}") for i in range(len(string_list))]
    return string_list

f = open("../_pdfs_2close/DR144.4_short.asc","r")
string = f.read()
zlines = string.split('\n')   # Split in to \n
# def read_info(ID):
    # line = [zlines[i] for i in range(len(zlines)) if ID in zlines[i]]
    # if line != []:
    #     z = float(line[0].split(' ')[-1])
    # else:
    #     z = -99
    # return z

from astropy.cosmology import FlatLambdaCDM
# cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

print("ID & RA & DEC & z & Separ.& Mag.& Grism & Lines & PA\\\\")
print("&&&&($''$, kpc)&(pair)&&&(deg)\\\\")
offset_kpc_h0_list, z_list = [], []
for ID in ID_list:
    show_ID = ID[:4] + ID[9:14]
    line = [zlines[i] for i in range(len(zlines)) if ID in zlines[i]]
    line[0] = line[0].replace('  ', ' ')
    z = float(line[0].split(' ')[-1])
    RA, Dec = line[0].split(' ')[1], line[0].split(' ')[2]
    
    files_1 = glob.glob('../proof2close_HSC_images_5band/*/' + ID + '/fit_result/')
    files_2 = glob.glob('../extra/*/fit_result*/' + ID + '/')
    files = files_1 + files_2
    file = glob.glob(files[-1]+'fit_result_{0}-band.txt'.format('I'))
    if file != []:
        f = open(file[0],"r")    
    string = f.read()    
    lines = string.split('\n')   # Split in to \n
    trust = 2    
    l1 = [i for i in range(len(lines)) if 'model_PS_result:' in lines[i]]
    AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])    
    AGN_pos = np.array([[-1*AGN_dic[i]['ra_image'], AGN_dic[i]['dec_image']] for i in range(len(AGN_dic))])    
    offset = np.sum( (AGN_pos[0] -  AGN_pos[1])**2)**0.5 
    scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
    offset_kpc = offset * scale_relation   #In kpc
    Mags = [AGN_dic[0]['magnitude'], AGN_dic[1]['magnitude']]
    print(show_ID, RA, Dec, z, '{0:.2f}, {1:.2f}'.format(offset, offset_kpc), 
          '{0:.2f}, {1:.2f}'.format(np.min(Mags), np.max(Mags))  
          )
    offset_kpc_h0_list.append(offset_kpc * 70 / 100)
    z_list.append(z)
    
#%%
import sys
sys.path.insert(0,'Shenli_materials/shenli_figure1/')
import separation
plt.scatter(np.array(offset_kpc_h0_list), np.array(z_list),c = 'red', marker = '*', s = 25, alpha = 0.9, label = 'Our sample')
plt.legend(loc='upper left', prop={'size': 7},ncol=1)
plt.show()