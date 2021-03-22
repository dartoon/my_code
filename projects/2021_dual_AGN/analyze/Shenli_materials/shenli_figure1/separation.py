#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 21:50:25 2020

@author: tang
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pandas as pd
import os
from Ned_calculator import Cosmic

data_path = '/Users/Dartoon/Astro/Projects/my_code/projects/2021_dual_AGN/analyze/Shenli_materials/shenli_figure1/'
os.chdir(data_path)
# sample = pd.read_csv('whole_sample.csv', index_col = 0)
# binaries = sample[(sample['status']=='QSO')]
#binaries_gemini = sample[(sample['status']=='QSO')&(sample['telescope']=='gemini')]
#print(binaries)

# convert angular separation to physical separation in kpc/h
def ang_to_d(ang, z):
    cosmic = Cosmic(z=z)
    return ang * 4.848137 * 1e-3 * cosmic.DA_Mpc * cosmic.Ho / 100
# print (ang_to_d(1, 1.5))
#%%
Ef_data = pd.DataFrame(np.loadtxt(data_path+'Efte.txt'), columns = ['Sep(")', 'zspec','v', 'R'])
He2010_data = pd.DataFrame(np.loadtxt(data_path+'Hennawi2010.txt'), columns = ['zspec', 'v', 'Sep(")', 'R'])
He2006_data = pd.DataFrame(np.loadtxt(data_path+'Hennawi2006.txt'), columns = ['Sep(")', 'zspec', 'v', 'R'])
#seren = pd.read_csv(data_path+'serendipitous.csv',header=0)
Kayo_data = pd.DataFrame(np.loadtxt(data_path+'Kayo.txt'), columns = ['z', 'Sep(")'])
More_data = pd.DataFrame(np.loadtxt(data_path+'More2016.txt'), columns = ['z', 'Sep(")'])
Kayo_data['R'] = Kayo_data.apply(lambda x: ang_to_d(x['Sep(")'], x['z']), axis=1)
More_data['R'] = More_data.apply(lambda x: ang_to_d(x['Sep(")'], x['z']), axis=1)
# binaries['R'] = binaries.apply(lambda x: ang_to_d(x['Sep(")'], x['Redshift']), axis=1)
#binaries_gemini['R'] = binaries_gemini.apply(lambda x: ang_to_d(x['Sep(")'], x['Redshift']), axis=1)
# keck = binaries[binaries['telescope']=='keck']
#print(promising['Redshift'])
#print(He2010_data['R'])
#%%
plt.rcParams['savefig.dpi'] = 300 #图片像素
plt.rcParams['figure.dpi'] = 150 #分辨率
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

plt.scatter(He2006_data['R'], He2006_data['zspec'],c = 'blue', marker = 'v', s = 5, alpha = 0.8, label = 'Hennawi+06')
plt.scatter(He2010_data['R'], He2010_data['zspec'],c = 'green', marker = '^', s = 5, alpha = 0.8, label = 'Hennawi+10')
#plt.scatter(seren['Sep(kpc)']/0.7, seren['Redshift'],c = 'grey', marker = 's', s = 5, alpha = 0.8, label = 'serendipitous')
plt.scatter(Kayo_data['R'], Kayo_data['z'],c = 'grey', marker = 's', s = 5, alpha = 0.8, label = 'Kayo+12')
plt.scatter(More_data['R'], More_data['z'],c = 'deeppink', marker = 'x', s = 10, alpha = 0.8, label = 'More+16')
plt.scatter(Ef_data['R'], Ef_data['zspec'],c = 'purple', marker = 'o', s = 5, alpha = 0.8, label = 'Eftekharzadeh+17')
#plt.scatter(keck['R'], keck['Redshift'], c = 'orange', marker = '*', edgecolors='black',s=60,linewidths=0.5, label = 'Silverman+20')
# plt.scatter(binaries['R'], binaries['Redshift'], c = 'lime', marker = '*', edgecolors='black',s=80,linewidths=0.6, label = 'This work')
#plt.scatter(binaries['R'], binaries['Redshift'], c = 'gold', marker = '*', edgecolors='black',s=80,linewidths=0.6, label = 'gemini sources')

plt.xscale("log")
plt.ylim(bottom = 0)
xmin, xmax, ymin, ymax = plt.axis()
# print(ymin,ymax)
z_array = np.linspace(ymin,ymax,num = 100)
R_array1 = np.array([ang_to_d(0.6, x) for x in z_array])
R_array2 = np.array([ang_to_d(3, x) for x in z_array])
#print(z_array,R_array1)

plt.plot(R_array1, z_array, 'black', linestyle = '--', linewidth = 0.7)
plt.plot(R_array2, z_array, 'black', linestyle = '--', linewidth = 0.7)
#plt.axvline(10,ls='--',lw=0.7,c='black')
plt.annotate('\u03B8 = 0.6"', xy=(3.5,4), size = 7)
plt.annotate('\u03B8 = 3"', xy=(20,4), size = 7)
# plt.legend(loc='upper left', prop={'size': 8})
plt.xlim(xmin=0.2)
plt.xlabel('R (kpc/h)',  size = 12)
plt.ylabel('z (redshift)', size = 12)
plt.tick_params(axis="y",direction="in", right=True)
plt.tick_params(axis="x",which = 'both', direction="in", top = True)
# os.chdir('../plots')
# plt.savefig('separation.png', dpi = 300, bbox_inches='tight')
# plt.show()

# for data in [Ef_data, He2010_data, He2006_data, Kayo_data]:
#     print(len(data['R']))

