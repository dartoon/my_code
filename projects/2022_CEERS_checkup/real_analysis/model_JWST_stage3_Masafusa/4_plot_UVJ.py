#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:23:12 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
# ID, mags, z = 'idx0', 
# 1,2,0,51,35
# folder = '202209' #Consider only JWST
# idx = [1,2,0,51,35]

# #0.2 mag error
# folder = '202209' #Not HST

# #0.4 mag error
# folder = '20220901_' #Not HST
# folder = '20220901' #HST upper limit
folder = '20220906' #mag error based on PSF lib

from functions_for_result import esti_smass, load_prop, load_info, name_list

color = []
target_id_list = []
for idx in [1,2,0,51,35]:
    target_id, z = load_info(idx)
    steller_file = glob.glob('esti_smass/'+folder+str(idx)+'/gsf_spec_*.fits')[0]
    hdul = pyfits.open(steller_file)
    info = hdul[1].header 
    uv = info['UV50']
    vj = info['VJ50']
    color.append([vj, uv])
    target_id_list.append(name_list[idx])

color = np.array(color)

plt.figure(figsize=(11,11))

empty = np.array([True, False, True, False, False])
confrim = empty==False
# empty = np.array([True, False, False, False, False])
plt.scatter(color[:,0][empty], color[:,1][empty],s=880,marker="H",
            # c='darkred',s=280,marker="o",zorder=1, vmin=0.5, vmax=5, edgecolors='white')
             facecolors='white', edgecolors='black', linewidths=2, zorder=4, alpha =0.8)
jwst_p = plt.scatter(color[:,0][confrim], color[:,1][confrim],
            c='lightskyblue',s=880,marker="H",zorder=2, vmin=0.5, vmax=5, edgecolors='black')

for i in range(len(target_id_list)):
    if target_id_list[i] == 'SDSS1419':
        plt.text(color[i,0]-0.6, color[i,1]+0.1,target_id_list[i],  fontsize = 21, zorder =2,
                  bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 3})   
    if target_id_list[i] != 'SDSS1419':
        plt.text(color[i,0]+0.15, color[i,1],target_id_list[i],  fontsize = 21, zorder =2,
                 bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 3})   
    
#plt.title('', fontsize=27)
#%%
import pickle
# pickle.dump(results, open('size_mass_CANDELS_result.pkl', 'wb'))
results = pickle.load(open('size_mass_CANDELS_result.pkl','rb'))

relation = 0  # 0 M*- Reff ; 1 M*-color; 2 M*-n; 3 Reff-n relation
# z_range = [1.2,1.7]
z_range = [1.6,3.5]
z_cut = ((results[:,0]>z_range[0]) * (results[:,0]<z_range[1]))

all_color = []
for i in range(len(results)):
    all_color.append([results[:,9][i][0], results[:,9][i][1], results[:,9][i][2]])
all_color = np.asarray(all_color)


ColorUV = -(2.5* np.log10(all_color[:,0])-2.5* np.log10(all_color[:,1]))
ColorVJ = -(2.5* np.log10(all_color[:,1])-2.5* np.log10(all_color[:,2]))

blue_galaxy, red_galaxy = [], []

k = 1.17
b = 0.3
for i in range(len(all_color)):
    if ColorVJ[i] < 0.8425:
        blue_galaxy.append((ColorUV[i] < 1.286))
        red_galaxy.append((ColorUV[i] > 1.286))
    else:
        line_p = k*ColorVJ[i]+b
        blue_galaxy.append((ColorUV[i] -line_p < 0))        
        red_galaxy.append((ColorUV[i] -line_p > 0))
blue_galaxy = np.asarray(blue_galaxy)
red_galaxy = np.asarray(red_galaxy)
# cmap_r = matplotlib.cm.get_cmap('RdBu_r')
c1 = plt.scatter(ColorVJ[z_cut* blue_galaxy],ColorUV[z_cut * blue_galaxy],
            c='lightskyblue',s=280,marker=".",zorder=-90, alpha=0.6, edgecolors='white', #cmap=cmap_r, 
            label='CANDELS galaxy, star-forming')
c2 = plt.scatter(ColorVJ[z_cut* red_galaxy],ColorUV[z_cut * red_galaxy],
            c='darksalmon',s=280,marker=".",zorder=-90, alpha=0.6, edgecolors='white', #cmap=cmap_r, 
            label='CANDELS galaxy, quiescent')
plt.plot( np.linspace(-1,0.8425), np.linspace(-1,0.8425)*0 + 1.286,'k--')
plt.plot( np.linspace(0.8425,2), k*np.linspace(0.8425,2)+b,'k--')
# plt.plot( np.linspace(-1,0.8425), k*np.linspace(-1,0.8425)+b)

#%%
plt.xlabel("V$-$J",fontsize=27)
plt.ylabel("U$-$V",fontsize=27)
plt.tick_params(labelsize=20)
plt.xlim([-1,3])
plt.ylim([0,2.5])
plt.tick_params(labelsize=25)
# from matplotlib.legend_handler import HandlerTuple
plt.legend([c1, c2, jwst_p], ['CANDELS galaxy, star-forming', 'CANDELS galaxy, quiescent', 
                              'This work, 1.6<z<3.5'],
               loc='upper left',fontsize=21,numpoints=1)

#plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.savefig('outcomes/UVJ.png')
plt.show()
# # 
# steller_file = glob.glob('esti_smass/'+folder+str(idx)+'/gsf_spec_*.fits')[0]
# hdul = pyfits.open(steller_file)
# info = hdul[0].header