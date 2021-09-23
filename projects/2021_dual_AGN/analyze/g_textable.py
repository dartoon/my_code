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

from tools import read_string_list, read_info, cal_oreination

f = open("material/ID_RA_DEC_z.txt","r")
string = f.read()
zlines = string.split('\n')   # Split in to \n
# def read_info(ID):
    # line = [zlines[i] for i in range(len(zlines)) if ID in zlines[i]]
    # if line != []:
    #     z = float(line[0].split(' ')[-1])
    # else:
    #     z = -99
    # return z
# Halpha = 6562.8

em_lines = np.array([1549, 2798, 4861, 5007, 6563])  #CIV, MgII, Hb, OIII, Ha
# G102range = [8000, 11500]
# G141range = [10750, 17000]
# def av_filter(z):
#     lines = np.array([CIV, MgII, Hb, OIII])
#     redshift_lines = (1+z) * lines
#     G102_bool =  (redshift_lines>G102range[0]+100) * (redshift_lines<G102range[1]-100)
#     G141_bool =  (redshift_lines>G141range[0]+100) * (redshift_lines<G141range[1]-100)
#     # return G102_bool, G141_bool
#     s1 = np.array(['CIV', 'MgII', 'H$beta$', '[OIII]'])[G102_bool] 
#     s2 = np.array(['CIV', 'MgII', 'H$beta$', '[OIII]'])[G141_bool] 
#     s1 = [s1[i] for i in range(len(s1))]
#     s2 = [s2[i] for i in range(len(s2))]
#     # str1 = "G102: " + repr(s1)
#     # str2 = " G141: " + repr(s2)    
#     # s = str1 + str2
#     if s2 != []:
#         try:
#             s = "G141 & " + s2[0] + '+' + s2[1]
#         except:
#             s = "G141 & " + s2[0]
#     elif s2 == [] and s1 != []:
#         try: 
#             s = "G102 & " + s1[0] + '+' +  s1[1]
#         except:
#             s = "G102 & " + s1[0]
#     else:
#         s = "No fileter!!! & "
#     return s

from astropy.cosmology import FlatLambdaCDM
# cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

print("ID & RA & DEC & z & Separ.& Mag.& PA & lines (CIV, MgII, Hb, OIII, Ha )\\\\")
print("&&&&($''$, kpc)&(pair)&&&(deg)\\\\")
offset_kpc_h0_list, z_list = [], []
for ID in ID_list:
    show_ID = ID[:4] + '$' + ID[9] + '$' + ID[10:14]
    line = [zlines[i] for i in range(len(zlines)) if ID in zlines[i]]
    line[0] = line[0].replace('  ', ' ')
    z = float(line[0].split(' ')[-1])
    RA, Dec = line[0].split(' ')[1], line[0].split(' ')[2]
    # RA, Dec, z = read_info(ID)
    
    files_1 = glob.glob('../proof2close_HSC_images_5band/*/' + ID + '/fit_result/')
    files_2 = glob.glob('../extra/*/fit_result*/' + ID + '/')
    files = files_1 + files_2
    file = glob.glob(files[-1]+'fit_result_{0}-band.txt'.format('I'))
    if file != []:
        f = open(file[0],"r")    
    string = f.read()    
    lines = string.split('\n')   # Split in to \n
    trust = 2    
    obs_lines = em_lines * (1+z)/ 1.e4 #Units of um
    l1 = [i for i in range(len(lines)) if 'model_PS_result:' in lines[i]]
    AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])    
    AGN_pos = np.array([[-1*AGN_dic[i]['ra_image'], AGN_dic[i]['dec_image']] for i in range(len(AGN_dic))])    
    offset = np.sum( (AGN_pos[0] -  AGN_pos[1])**2)**0.5 
    scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
    offset_kpc = offset * scale_relation   #In kpc
    Mags = [AGN_dic[0]['magnitude'], AGN_dic[1]['magnitude']]
    print(show_ID, '&' , RA, '&' , Dec, '&' , '{0:.3f}'.format(z), '&' , '{0:.2f},{1:.1f}'.format(offset, offset_kpc), '&', 
          '{0:.1f},{1:.1f}'.format(np.min(Mags), np.max(Mags)), '& {0:.1f}'.format(cal_oreination(ID)), '& {0} $  \\\\'.format( np.around(obs_lines,3) ) )
    offset_kpc_h0_list.append(offset_kpc * 70 / 100)
    z_list.append(z)
    
# #%%
# import sys
# sys.path.insert(0,'Shenli_materials/shenli_figure1/')
# import separation
# import matplotlib as mat
# mat.rcParams['font.family'] = 'STIXGeneral'

# plt.scatter(np.array(z_list), np.array(offset_kpc_h0_list), c = 'red', marker = '*', edgecolors='black', s = 305, alpha = 0.9, label = 'Proposed Sample')

# import pandas as pd
# shenli_sample = pd.read_csv('../whole_sample_new.csv', index_col = 0)
# shenli_sep = shenli_sample['Sep(")']
# shenli_z = shenli_sample['Redshift']
# shenli_tel = shenli_sample['telescope']
# shenli_stat = shenli_sample['status']
# shenli_sep_l, shenli_z_l = [], []
# for i in range(len(shenli_sample)):
#     if shenli_stat[i] == 'QSO': 
#         scale_relation = cosmo.angular_diameter_distance(shenli_z[i]).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
#         shenli_sep_l.append(shenli_sep[i]* scale_relation )
#         shenli_z_l.append(shenli_z[i])
# shenli_sep_l = np.array(shenli_sep_l)
# shenli_z_l = np.array(shenli_z_l)


# plt.scatter(0.2, 0.430* 70 / 100,c = 'blue', marker = 'o', edgecolors='black', s = 50, alpha = 0.9, label = 'Goulding+19')

# plt.scatter(shenli_z_l, shenli_sep_l * 70 / 100, marker="h",edgecolors='black',
#             c='m', s=220,zorder=10,alpha=1, label = 'Our paper')


# DeRosa = np.array([[0.0749, 30], [0.0551, 43], [0.0482, 51], [0.0446, 59] ]) #De Rosa MNRAS 2018
# plt.scatter(DeRosa[:,0], DeRosa[:,1] * 70 / 100, marker="v",  edgecolors='black',
#     c='green', s=55,zorder=10,alpha=1, label = 'De Rosa+18')

# plt.scatter(0.858, 2 * 70 / 100, marker="^",  edgecolors='black',
#     c='blue', s=55,zorder=10,alpha=1, label = 'Shields+2012')

# plt.plot(np.linspace(0,1000)*0+1,np.linspace(0,1000), '--', c = 'black', linewidth = 1.5)

# plt.text(0.2, 400, r"Confirmed dual QSOs in the literature",fontsize=25, color='black', bbox = {'facecolor':'white','alpha': 0.5} )

# # 233713.66+005610.8: pos1 = np.array([0.139, -0.064]) pos2 = np.array([0.474, 1.227]) np.sum((pos1 - pos2)**2)**0.5 
# handles, labels = plt.gca().get_legend_handles_labels()
# order = [5,7,0,1,2,3,6,8,4,9]
# plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc=4, prop={'size': 16}, ncol=2, frameon=True,  
#            bbox_to_anchor=(1, -0.01))
# plt.tick_params(labelsize=20)
# ax = plt.axes()
# ax.set_xticks(np.arange(0, 5, 0.5))

# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)#, color='r’)
# # plt.savefig('offset_vs_z.png')
# plt.show()
# # mv offset_vs_z.png ../../../../../../../../Astro/proposal/2021_HST_Grism/