#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 21:26:34 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle

# sim_files = glob.glob("sim_result/*_QuickResult.pkl")
sim_files = glob.glob("sim_result_bulge_n2/*_QuickResult.pkl")

B2T_result_list, chisq_list, n_list, bulge_re_list, disk_re_list, flux_single_list = [], [], [], [],[], []
for i in range(len(sim_files)):
    result = pickle.load(open(sim_files[i],'rb'))
    bic_result, chisq_result, B2T_result, True_param, infer_param_single, infer_param_DB = result
    # if B2T_result[2] < 0.95 and chisq_result[-1]<1.4:
    if B2T_result[2] < 0.90 and chisq_result[-1]<3:# and infer_param_DB[1]['R_sersic'] - True_param[1]['R_sersic'] < 0.4 :
    # if 2>1:
        B2T_result_list.append(B2T_result[1:])
        chisq_list.append(chisq_result[1:])
        n_list.append([True_param[1]['n_sersic'], infer_param_DB[1]['n_sersic']])
        bulge_re_list.append([True_param[1]['R_sersic'], infer_param_DB[1]['R_sersic']])
        disk_re_list.append([True_param[2]['R_sersic'], infer_param_DB[2]['R_sersic']])
        flux_single_list.append([True_param[1]['flux_within_frame'] + True_param[2]['flux_within_frame'] , 
                                 infer_param_single[1]['flux_within_frame'],
                                 infer_param_DB[1]['flux_within_frame'] + infer_param_DB[2]['flux_within_frame']
                                 ])
    
B2T_results = np.array(B2T_result_list)
chisqs = np.array(chisq_list)
n_list = np.array(n_list)
bulge_re_list = np.array(bulge_re_list)
disk_re_list = np.array(disk_re_list)
flux_single_list = np.array(flux_single_list)

plt.figure(figsize=(9, 6))
plt.scatter(chisqs[:,1], (B2T_results[:,1]- B2T_results[:,0]) * 100 )
plt.xlabel("Fitted Reduced Chisq",fontsize=27)
plt.ylabel("B/T % (infer - truth)",fontsize=27)
plt.tick_params(labelsize=20)
plt.title(sim_files[0].split('/')[0][-8:], fontsize=27)
plt.xlim(0.8, 5)
plt.ylim(-40, 40)
plt.show()

#%%

plt.figure(figsize=(6, 5))
plt.title(sim_files[0].split('/')[0][-8:], fontsize=27)
plt.hist((B2T_results[:,1]- B2T_results[:,0])*100 )
plt.xlabel("B/T % (infer - truth)",fontsize=27)
plt.tick_params(labelsize=20)
plt.show()

plt.figure(figsize=(9, 6))
plt.scatter(chisqs[:,1], (n_list[:,1]- n_list[:,0]) )
plt.xlabel("Fitted Reduced Chisq",fontsize=27)
plt.ylabel("Sersic_n (infer - truth)",fontsize=27)
plt.tick_params(labelsize=20)
plt.title(sim_files[0].split('/')[0][-8:], fontsize=27)
plt.xlim(0.8, 5)
# plt.ylim(-4, 4)
plt.show()

plt.figure(figsize=(9, 6))
plt.scatter(chisqs[:,1], (bulge_re_list[:,1]- bulge_re_list[:,0]) )
plt.xlabel("Fitted Reduced Chisq",fontsize=27)
plt.ylabel("Bulge_reff (infer - truth)",fontsize=27)
plt.tick_params(labelsize=20)
plt.title(sim_files[0].split('/')[0][-8:], fontsize=27)
plt.xlim(0.8, 5)
# plt.ylim(-4, 4)
plt.show()

# plt.figure(figsize=(9, 6))
# plt.scatter(chisqs[:,1], (flux_single_list[:,2]- flux_single_list[:,0]) )
# plt.xlabel("Fitted Reduced Chisq",fontsize=27)
# plt.ylabel("infer total flux (infer - truth)",fontsize=27)
# plt.tick_params(labelsize=20)
# plt.title(sim_files[0].split('/')[0][-8:], fontsize=27)
# plt.xlim(0.8, 5)
# # plt.ylim(-4, 4)
# plt.show()



# plt.hist(4- n_list[:,1] )
# plt.scatter(chisqs[:,1], (4- n_list[:,1] ) )

# plt.scatter(bulge_re_list[:,0]/disk_re_list[:,0], B2T_results[:,0]- B2T_results[:,1] )
