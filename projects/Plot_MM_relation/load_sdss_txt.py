#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 16:05:10 2020

@author: dartoon
"""

import numpy as np
line_means = ['id', 'z', 'ra', 'dec', 'fix_sersic_n', 'sersic_n_fitted', 'sersic_re_fitted', 'sersic_n_corrected',
         'sersic_re_corrected', 'host_mag_g', 'host_mag_r', 'host_mag_i', 'host_mag_z', 'host_mag_y',
         'ps_mag_g', 'ps_mag_r', 'ps_mag_i', 'ps_mag_z', 'ps_mag_y', 'decomposition_chisq', 'stellar_mass', 
         'sed_chisq', 'logMBH', 'logMBH_err']
infers  = np.loadtxt('./sdss_quasar_decomposition_v1.txt', dtype=str)
IDs_ = infers[:, 0]
HSC_z_ = infers[:,1].astype(np.float)
HSC_Mstar_ = infers[:,20].astype(np.float)
HSC_MBHs_ = infers[:,22].astype(np.float)
HSC_MBHs_err_ = infers[:,23].astype(np.float)

HSC_RA_ = infers[:,2].astype(np.float)
HSC_DEC_ = infers[:,3].astype(np.float)

flags_  = np.loadtxt('./sdss_quasar_decomposition_v1_catalog_flag.txt', dtype=str)
flags = flags_[:,0]

IDs, HSC_z, HSC_Mstar, HSC_MBHs, HSC_MBHs_err =[], [], [], [], []
HSC_RA, HSC_DEC  = [], []
for i in range(len(IDs_)):
    idx = np.where(IDs_[i] == flags)[0][0]
    if flags_[idx][1] == 'y':
        IDs.append(IDs_[i])
        HSC_z.append(HSC_z_[i])
        HSC_Mstar.append(HSC_Mstar_[i])  #Uncertainty as 0.2 dex
        HSC_MBHs.append(HSC_MBHs_[i])
        HSC_MBHs_err.append(HSC_MBHs_err_[i])
        HSC_RA.append(HSC_RA_[i])
        HSC_DEC.append(HSC_DEC_[i])
    #%%
# import matplotlib.pyplot as plt  
# plt.figure(figsize=(11.5,12))      
# plt.scatter(HSC_Mstar,HSC_MBHs,c='blue',
#             s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.4)

# plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
# plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
# plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
# plt.xlim(9,12.5)
# plt.ylim(6.0,10.3)
# plt.grid(linestyle='--')
# plt.tick_params(labelsize=25)
# plt.show()
#Write files:
write_file = open('sdss_quasar_z_0.2-0.3.txt','w') 
write_file.write("#id ra dec z\n")
for i in range(len(HSC_z)):
    if HSC_z[i]>0.2 and HSC_z[i]<0.3:
        write_file.write("{0} {1} {2} {3}\n".format(IDs[i], HSC_RA[i], HSC_DEC[i], HSC_z[i]))
write_file.close()
        