#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 10:34:23 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from astropy.coordinates import SkyCoord
from astropy import units as u
from tools import read_string_list, read_z, cal_offset

# pick = True
# if pick == True:
# folders_1 = glob.glob('../proof2close_HSC_images_5band/*/'+ '*/')
# IDs_1 = [folders_1[i].split('/')[-2] for i in range(len(folders_1))]
# folders_2 = glob.glob('../proofBHBH/model_Iband_zover1/' + '*/')
# IDs_2 = [folders_2[i].split('/')[-2] for i in range(len(folders_2))]
# IDs = IDs_1 #+ IDs_2
# IDs = list(dict.fromkeys(IDs))
z_thre = 0
# else:
#     folders = glob.glob('../_John_fitted/*')  #All John's detected source.
#     IDs = [folders[i].split('/')[-1].split('_HSC')[0] for i in range(len(folders))]
#     z_thre = 0

from ID_list import ID_list as run_ID_list
offset_list, z_list, ID_list, mags_list = [], [], [], []

for ID in run_ID_list:
    string = ID
    folder_1 = glob.glob('../proof2close_HSC_images_5band/*/'+ ID+ '/')
    if folder_1 != []: # and 'z_below1' not in folder_1[0]:
        folder = folder_1[0] + 'fit_result/'
        file = folder + 'fit_result_I-band.txt'
    # elif folder_1 != [] and 'z_below1' in folder_1[0]:
    #     folder = '../_John_fitted/'+ID+'_HSC-I/'   #For these z below 1(or z unkonwn), not fitted and use John's fit.
    #     file = folder + 'fit_result.txt'
    else:
        folder_2 = glob.glob('../proofBHBH/model_Iband_zover1/'+ ID + '*/') 
        folder = folder_2[0]
        file = folder + 'fit_result_I-band.txt'
    f = open(file,"r")    
    string = f.read()
    lines = string.split('\n')   # Split in to \n
    l0 = [i for i in range(len(lines)) if 'PS PS center offset' in lines[i]]
    l1 = [i for i in range(len(lines)) if 'AGN mag:' in lines[i]]
    if l0 != []:
        offset = float(lines[l0[-1]].split(' ')[-1])
        z = read_z(ID)
        mags = lines[l1[-1]].split(' ')[2:4]
        mags = [float(mags[i]) for i in range(2)]
        if z > z_thre: #and offset>0.4:
            offset_list.append(offset)
            z_list.append(z)
            ID_list.append(ID)
            mags_list.append(mags)

#%%    
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
offset_list = np.array(offset_list)
z_list = np.array(z_list)
ID_list = np.array(ID_list)
mags_list = np.array(mags_list)

plt.figure(figsize=(11, 9))

p_ID_list = []
legend = 'Proposed sample'
# if pick == True:
ct = 0
special_ID = []
for ID in run_ID_list:
    k = [i for i in range(len(ID_list)) if ID == ID_list[i]]
    if k != []:
        k = k[0]
        plt.scatter(z_list[k], offset_list[k],
                c='red',s=405,marker="*",zorder=100, edgecolors='black',alpha=0.8, label = legend)
        legend = None
        ct = ct + 1
        p_ID_list.append(ID)
        special_ID.append(ID)
            
import pandas as pd
shenli_sample = pd.read_csv('Shenli_materials/whole_sample_new.csv', index_col = 0)
shenli_ID = shenli_sample['ID']
shenli_sep = shenli_sample['Sep(")']
shenli_qsog = shenli_sample['qso_g']
shenli_qsor = shenli_sample['qso_r']
shenli_comg = shenli_sample['com_g']
shenli_comr = shenli_sample['com_r']
shenli_z = shenli_sample['Redshift']
shenli_tel = shenli_sample['telescope']
shenli_stat = shenli_sample['status']
legends = ['Confirmed QSO pairs', 'Confirmed nonQSO pairs', 'Gemini/NTT confirming', 'xxx']
for i in range(len(shenli_sample)):
    offset = cal_offset(shenli_ID[i])
    if offset == None:
        print(shenli_ID[i], "use Shenli")
        offset = shenli_sep[i]
    if shenli_stat[i] == 'QSO': 
        plt.scatter(shenli_z[i], offset, marker="h",edgecolors='black',
            c='m', s=220,zorder=10,alpha=1, label = legends[0])
        legends[0] = None
        p_ID_list.append(shenli_ID[i])
    labels = ['QSO', 'wait', 'done', 'issued']
    if shenli_stat[i] not in labels: 
        plt.scatter(shenli_z[i], offset, marker="X",
            c='green', edgecolors='black', s=150,zorder=10,alpha=1, label = legends[1])
        legends[1] = None
        p_ID_list.append(shenli_ID[i])
        if shenli_sep[i]<0.5:
            print(shenli_ID[i])
    if shenli_stat[i] == 'done' or shenli_stat[i] == 'issued': 
        plt.scatter(shenli_z[i], offset, marker="^",
            c='c', edgecolors='black',s=100,zorder=9,alpha=0.9, label = legends[2])
        legends[2] = None
        if shenli_ID[i] in run_ID_list:
            print(shenli_ID[i])
        p_ID_list.append(shenli_ID[i])
    # if shenli_stat[i] == 'wait': 
    #     if (shenli_comg[i] - shenli_comr[i] ) <1 and (shenli_qsog[i] - shenli_qsor[i] ) <1:
    #         plt.scatter(shenli_z[i], shenli_sep[i],
    #             c='black',s=200,marker=".",zorder=-10,alpha=0.6, label = legends[3])
    #         legends[3] = None

file_all_cand_0 = glob.glob('../proofBHBH/allband_data/fit_result/*/')
file_all_cand_1 = glob.glob('../proof2close_HSC_images_5band/z_over1/*/')
file_all_cand =  file_all_cand_0 + file_all_cand_1
legend = 'All Candidates'

for i in range(len(file_all_cand)):
# for i in range(2):    
    folder = file_all_cand[i]
    ID = file_all_cand[i].split('/')[-2]
    if 'proof2close_HSC_images_5band' in folder:
        folder = folder + 'fit_result' + '/'
    # else:
    #     ID = file_all_cand[i].split('/')[-2]
    if ID not in p_ID_list:
        trust = 2    
        if_p = True
        bands = ['G', 'R']
        AGN_mags = []
        for band in bands:
            file = glob.glob(folder + 'fit_result_{0}-band.txt'.format(band))
            if file != []:
                f = open(file[0],"r")    
                string = f.read()    
                lines = string.split('\n')   # Split in to \n
                l1 = [i for i in range(len(lines)) if 'model_PS_result:' in lines[i]]
                AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])    
                AGN_mags.append([AGN_dic[0]['magnitude'], AGN_dic[1]['magnitude']])
            else:
                if_p = False
        if if_p == True:
            if AGN_mags[0][0] - AGN_mags[1][0]<1 and AGN_mags[0][1] - AGN_mags[1][1]<1: 
                ID_z = read_z(ID)
                offset = cal_offset(ID)
                if ID_z<4:
                    plt.scatter(ID_z, offset, marker="o",edgecolors='black',
                        c='lightsalmon',s=80, zorder=-10,alpha=0.5, label = legend)
                    if ID_z>1. and offset<0.7:
                        special_ID.append(ID)
                # if ID_z>4:
                #     print(ID) #result: 022657.63-033335.4
                legend = None
                p_ID_list.append(ID)
                
# special_ID = special_ID+ run_ID_list                
# special_ID = list(dict.fromkeys(special_ID))
# len(special_ID) #result:39

plt.plot(np.linspace(-0.2,6)*0 + 1 ,np.linspace(-0.2,6), '--' , c ='black', linewidth = 2.0 )
plt.xlabel("Redshift",fontsize=27)
plt.ylabel("Projected separation (arcsec)",fontsize=27)
plt.tick_params(labelsize=20)
plt.plot(np.linspace(0,5),np.linspace(0,5)*0+0.7, '--' ,c = 'black', linewidth = 1.5   )
# plt.plot(np.linspace(0,1)*0 + 1.3 ,np.linspace(-0.2,0.7), '--', c = 'black', linewidth = 1.5   )
plt.plot(np.linspace(0,1)*0 + 2.35 ,np.linspace(-0.2,0.7), '--', c = 'black' , linewidth = 1.5  )
plt.plot(np.linspace(0,1)*0 + 3.1 ,np.linspace(-0.2,0.7), '--', c = 'black' , linewidth = 1.5  )
# plt.text(0.1, -0.08, r"G102 H$\beta$+[OIII]",fontsize=19, color='orange')
plt.text(1.3, -0.09, r"G141 H$\beta$+[OIII]",fontsize=19, color='blue')
plt.text(2.4, -0.09, r"G102 MgII",fontsize=19, color='blue')
plt.text(3.2, -0.09, r"G141 MgII",fontsize=19, color='blue')
plt.text(4.0, 0.4, r"HST",fontsize=29, color='green')
plt.text(3.3, 0.9, r"Ground-based telescope",fontsize=19, color='green')
plt.text(0.2, 4.4, r"HSC dual QSO Candidates",fontsize=25, color='black', bbox = {'facecolor':'white','alpha': 0.5} )
plt.ylim(-0.2, 5)
plt.xlim(0., 4.7)
handles, labels = plt.gca().get_legend_handles_labels()
order = [0,3,1,2]

ax = plt.axes()
plt.tick_params(axis="x",direction="in", right=True)
plt.tick_params(axis="y",which = 'both', direction="in", top = True)
plt.tick_params(which='both', width=1)
plt.tick_params(which='major', length=7)
ax.set_xticks(np.arange(0, 5, 0.5))
ax.set_yticks([0,0.7,1,2,3,4,5])
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper right', prop={'size': 20},ncol=1)
# plt.savefig('sample_select.png')
plt.show()

special_ID = list(dict.fromkeys(special_ID))

#%%
# plt_ID_list = ID_list[(np.max(mags_list,axis=1)<23)*(offset_list<1)] 
# for i in range( len(plt_ID_list ) ):
#     ID = plt_ID_list[i]
#     if ID not in run_ID_list:
#         print(ID)
    
    