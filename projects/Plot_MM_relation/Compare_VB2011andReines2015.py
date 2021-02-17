#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:02:45 2021

@author: Dartoon
"""
import numpy as np
import matplotlib.pyplot as plt

f = open('data/Reines_2015_table_1.txt',"r")
# ID  Running identification number
# NSAID NASA-Sloan Atlas identification number
# SDSS SDSS identification
# P-MJD-F Plate-MJD-Fiber of SDSS spectrum
# z  Redshift (1)
# mag  iMag Absolute i-band magnitude (2)
# mag  g-i  The (g-i) color (2)
# [solMass] logM* Log stellar mass; corrected for AGN contribution (3)
# [solMass] logMBH Log black hole mass (3)
string = f.read()
Reines_t1 = string.split('\n')   # Split in to \n
Reines_t1 = [Reines_t1[i].split(' ') for i in range(len(Reines_t1)) if Reines_t1[i][0]!= '#']
Reines_t2 = np.loadtxt('data/Reines_2015_table_2.txt')
#   1-  3 I3     ---          ID        Running identification number
#   5- 11 F7.1   10-17mW/m2   Hb-n      Narrow H{beta} emission line flux (1)
#  13- 18 F6.1   10-17mW/m2 e_Hb-n      Error on Hb-n (1)
#  20- 26 F7.1   10-17mW/m2   Hb-b      Broad H-beta emission line flux (1)
#  28- 34 F7.1   10-17mW/m2 e_Hb-b      Error on Hb-b (1)
#  36- 42 F7.1   10-17mW/m2   OIII-5007 The [OIII] 5007 emission line flux (1)
#  44- 49 F6.1   10-17mW/m2 e_OIII-5007 Error on OIII-5007 (1)
#  51- 57 F7.1   10-17mW/m2   NII-6548  The [NII] 6548 emission line flux (1)
#  59- 64 F6.1   10-17mW/m2 e_NII-6548  Error on NII-6548 (1)
#  66- 72 F7.1   10-17mW/m2   Ha-n      Narrow H{alpha} emission line flux (1)
#  74- 79 F6.1   10-17mW/m2 e_Ha-n      Error on Ha-n (1)
#  81- 88 F8.1   10-17mW/m2   Ha-b      Broad H{alpha} emission line flux (1)
#  90- 95 F6.1   10-17mW/m2 e_Ha-b      Error on Ha-b  (1)
#  97-103 F7.1   10-17mW/m2   NII-6583  The [NII] 6583 emission line flux (1)
# 105-110 F6.1   10-17mW/m2 e_NII-6583  Error on NII-6583 (1)
# 112-118 F7.1   10-17mW/m2   SII-6716  The [SII] 6716 emission line flux (1)
# 120-125 F6.1   10-17mW/m2 e_SII-6716  Error on SII-6716 (1)
# 127-133 F7.1   10-17mW/m2   SII-6731  The [SII] 6731 emission line flux (1)
# 135-140 F6.1   10-17mW/m2 e_SII-6731  Error on SII-6731 (1)
# 142-148 F7.1   km/s         FWHM      Full width at half maximum 
f = open('data/Reines_2015_table_3.txt',"r")
# 1 A1 --- Type Galaxy type (1)
# 3- 16 A14 --- Name Object name (2)
# 18 A1 --- f_Name Flag on Name (3)
# 20- 24 F5.2 [solMass] logM* Log total stellar mass (4)
# 26- 30 F5.2 [solMass] logMBH Log black hole mass (5)
# 33- 36 F4.2 [solMass] e_logMBH Error in logMBH
string = f.read()
Reines_t3 = string.split('\n')   # Split in to \n
Reines_t3 = [Reines_t3[i].split(' ') for i in range(len(Reines_t3)) if Reines_t3[i][0]!= '#']

#%% Load Vardha 2011
f1 ='data/Bennert+2011_local.txt'
b11_l = np.loadtxt(f1)[:,1:]
# b11_local_Reff = b11_l[:,8]
b11_local_MBH = b11_l[:,6]
b11_local_mstar_sph = b11_l[:,4]  #spheroid stellar mass
b11_local_mstar = b11_l[:,9]  #!!! Change to total mass

#%%
Reines_t1_Mstar = [float(Reines_t1[i][-2]) for i in range(len(Reines_t1))]
# Reines_t1_MBH = [float(Reines_t1[i][-1]) for i in range(len(Reines_t1))]
Reines_t1_MBH = [(float(Reines_t1[i][-1]) -0.2) for i in range(len(Reines_t1))]

Reines_t3_g1_Mstar = [float(Reines_t3[i][-3]) for i in range(len(Reines_t3)) if float(Reines_t3[i][0]) == 1]
Reines_t3_g1_MBH = [float(Reines_t3[i][-2]) for i in range(len(Reines_t3)) if float(Reines_t3[i][0]) == 1]
Reines_t3_g2_Mstar = [float(Reines_t3[i][-3]) for i in range(len(Reines_t3)) if float(Reines_t3[i][0]) == 2]
Reines_t3_g2_MBH = [float(Reines_t3[i][-2]) for i in range(len(Reines_t3)) if float(Reines_t3[i][0]) == 2]
Reines_t3_g3_Mstar = [float(Reines_t3[i][-3]) for i in range(len(Reines_t3)) if float(Reines_t3[i][0]) == 3]
Reines_t3_g3_MBH = [float(Reines_t3[i][-2]) for i in range(len(Reines_t3)) if float(Reines_t3[i][0]) == 3]

#%%
plt.figure(figsize=(11.5,12))
plt.scatter(b11_local_mstar,b11_local_MBH,s=180, c ='red',
            marker="o",zorder=100,  edgecolors='white', label='Bennert 2011, Mstellar total mass')     
plt.scatter(b11_local_mstar_sph,b11_local_MBH,s=180, c ='yellow',
            marker="o",zorder=99,  edgecolors='white', label='Bennert 2011, Mstellar spheroid mass')     
plt.scatter(Reines_t1_Mstar,Reines_t1_MBH, s=200, c ='blue',
            marker="o",zorder=95,  edgecolors='black', label='Reines 2015, board-line AGN')
plt.scatter(Reines_t3_g1_Mstar,Reines_t3_g1_MBH, s=200, c = 'limegreen',
            marker="s",zorder=95,  edgecolors='black', label='Reines 2015, board-line AGN table3')
plt.scatter(Reines_t3_g2_Mstar,Reines_t3_g2_MBH, s=200, c = 'violet',
            marker="s",zorder=95,  edgecolors='black', label='Reines 2015, RM AGN in table 3')
plt.scatter(Reines_t3_g3_Mstar,Reines_t3_g3_MBH, s=200, c = 'olive',
            marker="s",zorder=95,  edgecolors='black', label='Reines 2015, inactive galaxies in table3')


m_ml, b_ml = (0.981, -2.5459)
xl = np.linspace(5, 13, 100)
plt.plot(xl, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)

plt.xlim(8.5,12.5)
plt.ylim(4.0,11)
plt.tick_params(labelsize=25)
plt.legend(prop={'size':15}, loc = 2)
plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
plt.show()
