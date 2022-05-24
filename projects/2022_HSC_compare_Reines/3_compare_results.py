#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 17:30:14 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
Reines_t1 = np.loadtxt('./2021_previous/Reines_2015_table_1.txt', dtype=str)
IDs = ['J131310.12+051942.1','J233837.09-002810.4',
       'J145819.55+045451.7', 'J143450.63+033842.5', 'J141920.64+043623.3', 'J141630.82+013708.0',
       'J140018.41+050242.2', 'J134426.41+441620.0','J131305.81+012755.9','J121826.72-000750.1',
       'J104252.93+041441.1','J012159.81-010224.3','J084143.50+013149.8','J095540.47+050236.5']

masses = []
for ID in IDs:
    idx = np.where(Reines_t1[:,2] == ID)[0][0]
    Reines_iMag = float(Reines_t1[idx, 5])
    z = float(Reines_t1[idx, 4])
    Reines_smass = float((Reines_t1[Reines_t1[:,2]==ID] )[0][-2])
    
    filename_p  = 'esti_smass/'+ID+'/SFH_*.fits'
    filename_p = glob.glob(filename_p)[0]
    hdul = pyfits.open(filename_p)
    table = hdul[1].data
    name = hdul[1].columns
    smass_idx = [i for i in range(len(name)) if 'Mstel50' in str(name[i])][0]
    inf_smass = table[1][smass_idx] # 'smass:'
    print(ID, Reines_smass, inf_smass)
    masses.append([Reines_smass, inf_smass])
masses = np.array(masses)
plt.figure(figsize=(11, 11))
plt.scatter(masses[:,0],masses[:,1],
            c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
for i in range(len(IDs)):
    plt.text(masses[i,0],masses[i,1],IDs[i][:7],fontsize=15)
plt.xlabel("Reines stellar mass",fontsize=27)
plt.ylabel("HSC inferred stellar mass",fontsize=27)
plt.plot(np.linspace(9.5, 11.5), np.linspace(9.5, 11.5))
plt.tick_params(labelsize=20)
plt.xlim([9.5,11.5])
plt.ylim([9.5,11.5])
#plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.show()