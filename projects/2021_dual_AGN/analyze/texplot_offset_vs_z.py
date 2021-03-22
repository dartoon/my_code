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

pick = True

if pick == True:
    folders_1 = glob.glob('../proof2close_HSC_images_5band/*/'+ '*/')
    IDs_1 = [folders_1[i].split('/')[-2] for i in range(len(folders_1))]
    folders_2 = glob.glob('../proofBHBH/model_Iband/' + '*/')
    IDs_2 = [folders_2[i].split('/')[-2] for i in range(len(folders_2))]
    IDs = IDs_1 + IDs_2
    IDs = list(dict.fromkeys(IDs))
    z_thre = 0
else:
    folders = glob.glob('../_John_fitted/*')  #All John's detected source.
    IDs = [folders[i].split('/')[-1].split('_HSC')[0] for i in range(len(folders))]
    z_thre = 0

# f = open("../_pdfs_2close/DR144.4_short.asc","r")
f = open("material/ID_RA_DEC_z.txt","r")
string = f.read()
zlines = string.split('\n')   # Split in to \n
# write_file = open('ID_RA_DEC_z.txt','w') 

def read_z(ID):
    line = [zlines[i] for i in range(len(zlines)) if ID in zlines[i]]
    if line != []:
        # write_file.write(line[0] + '\n')
        z = float(line[0].split(' ')[-1])
    else:
        z = -99
    return z

offset_list, z_list, ID_list, mags_list = [], [], [], []
for ID in IDs:
    string = ID
    if pick == True:
        folder_1 = glob.glob('../proof2close_HSC_images_5band/*/'+ ID+ '/')
        if folder_1 != [] and 'z_below1' not in folder_1[0]:
            folder = folder_1[0] + 'fit_result/'
            file = folder + 'fit_result_I-band.txt'
        elif folder_1 != [] and 'z_below1' in folder_1[0]:
            folder = '../_John_fitted/'+ID+'_HSC-I/'   #For these z below 1(or z unkonwn), not fitted and use John's fit.
            file = folder + 'fit_result.txt'
        else:
            folder_2 = glob.glob('../proofBHBH/model_Iband/'+ ID + '*/') 
            folder = folder_2[0]
            file = folder + 'fit_result_I-band.txt'
    else:
        folder = '../_John_fitted/'+ID+'_HSC-I/'
        file = folder + 'fit_result.txt'  
    try:
        f = open(file,"r")    
    except:
        continue
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
    # z = read_z(ID)
    # if z>0:
    #     z_list.append(z)
# write_file.close()

#%%    
from ID_list import ID_list as run_ID_list
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

offset_list = np.array(offset_list)
z_list = np.array(z_list)
ID_list = np.array(ID_list)
mags_list = np.array(mags_list)

plt.figure(figsize=(11, 9))
# plt.scatter(z_list[np.max(mags_list,axis=1)<23], offset_list[np.max(mags_list,axis=1)<23],
#             c='black',s=280,marker=".",zorder=-10, edgecolors='white',alpha=0.7)
plt.scatter(z_list, offset_list,
            c='black',s=280,marker=".",zorder=-10, edgecolors='white',alpha=0.7)
if pick == True:
    ct = 0
    for ID in run_ID_list:
        k = [i for i in range(len(ID_list)) if ID == ID_list[i]]
        if k != []:
            k = k[0]
            plt.scatter(z_list[k], offset_list[k],
                    c='red',s=380,marker=".",zorder=0, edgecolors='white',alpha=1)
            # if offset_list[k]<0.5:
            #     print(ID, offset_list[k])
            ct = ct + 1
plt.xlabel("Redshift",fontsize=27)
plt.ylabel("Projected separation (arcsec)",fontsize=27)
plt.tick_params(labelsize=20)
plt.plot(np.linspace(0,5),np.linspace(0,5)*0+1, c = 'blue'  )
plt.plot(np.linspace(0,1)*0 + 1.3 ,np.linspace(-0.2,1), c = 'orange'  )
plt.plot(np.linspace(0,1)*0 + 2.35 ,np.linspace(-0.2,1), c = 'orange'  )
plt.plot(np.linspace(0,1)*0 + 3.1 ,np.linspace(-0.2,1), c = 'orange'  )
plt.text(0.6, 0.15, r"G102",fontsize=20, color='orange')
plt.text(0.6, -0.05, r"H$\beta$+[OIII]",fontsize=20, color='orange')
plt.text(1.5, 0.15, r"G141",fontsize=20, color='orange')
plt.text(1.5, -0.05, r"H$\beta$+[OIII]",fontsize=20, color='orange')
plt.text(2.4, 0.15, r"G102",fontsize=20, color='orange')
plt.text(2.4, -0.05, r"MgII",fontsize=20, color='orange')
plt.text(3.2, 0.15, r"G141",fontsize=20, color='orange')
plt.text(3.2, -0.05, r"MgII",fontsize=20, color='orange')
plt.ylim(-0.2, 5)
plt.xlim(0., 4.7)
# plt.savefig('show_material/sample_select.png')
plt.show()

#%%
# plt_ID_list = ID_list[(np.max(mags_list,axis=1)<23)*(offset_list<1)] 
# for i in range( len(plt_ID_list ) ):
#     ID = plt_ID_list[i]
#     if ID not in run_ID_list:
#         print(ID)
    
    