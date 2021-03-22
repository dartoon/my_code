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
    
CIV, MgII, Hb, OIII = 1549, 2798, 4861, 5007
G102range = [8000, 11500]
G141range = [10750, 17000]
def av_filter(z):
    lines = np.array([CIV, MgII, Hb, OIII])
    redshift_lines = (1+z) * lines
    G102_bool =  (redshift_lines>G102range[0]+100) * (redshift_lines<G102range[1]-100)
    G141_bool =  (redshift_lines>G141range[0]+100) * (redshift_lines<G141range[1]-100)
    # return G102_bool, G141_bool
    s1 = np.array(['CIV', 'MgII', 'H$beta$', '[OIII]'])[G102_bool] 
    s2 = np.array(['CIV', 'MgII', 'H$beta$', '[OIII]'])[G141_bool] 
    s1 = [s1[i] for i in range(len(s1))]
    s2 = [s2[i] for i in range(len(s2))]
    # str1 = "G102: " + repr(s1)
    # str2 = " G141: " + repr(s2)    
    # s = str1 + str2
    if s2 != []:
        try:
            s = "G141 & " + s2[0] + '+' + s2[1]
        except:
            s = "G141 & " + s2[0]
    elif s2 == [] and s1 != []:
        try: 
            s = "G102 & " + s1[0] + '+' +  s1[1]
        except:
            s = "G102 & " + s1[0]
    else:
        s = "No fileter!!! & "
    return s

from astropy.cosmology import FlatLambdaCDM
# cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

print("ID & RA & DEC & z & Separ.& Mag.& Grism & Lines & PA\\\\")
print("&&&&($''$, kpc)&(pair)&&&(deg)\\\\")
offset_kpc_h0_list, z_list = [], []
for ID in ID_list:
    show_ID = ID[:4] + '$' + ID[9] + '$' + ID[10:14]
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
    print(show_ID, '&' , RA, '&' , Dec, '&' , '{0:.3f}'.format(z), '&' , '{0:.2f},{1:.1f}'.format(offset, offset_kpc), '&', 
          '{0:.1f},{1:.1f}'.format(np.min(Mags), np.max(Mags)), '&', 
          av_filter(z), '& TBD \\\\' )
    offset_kpc_h0_list.append(offset_kpc * 70 / 100)
    z_list.append(z)
    
#%%
import sys
sys.path.insert(0,'Shenli_materials/shenli_figure1/')
import separation
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

plt.scatter(np.array(offset_kpc_h0_list), np.array(z_list),c = 'red', marker = '*', s = 25, alpha = 0.9, label = 'Our sample')
plt.legend(loc='upper left', prop={'size': 7},ncol=1)
# plt.savefig('offset_vs_z.png')
plt.show()
# mv offset_vs_z.png ../../../../../../../../Astro/proposal/2021_HST_Grism/