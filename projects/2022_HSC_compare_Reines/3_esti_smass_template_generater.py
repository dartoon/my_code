#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 16:21:29 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob, pickle
# import os

# Reines_t1 = np.loadtxt('./2021_previous/Reines_2015_table_1.txt', dtype=str)
# IDs = ['J131310.12+051942.1','J233837.09-002810.4',
#        'J145819.55+045451.7', 'J143450.63+033842.5', 'J141920.64+043623.3', 'J141630.82+013708.0',
#        'J140018.41+050242.2', 'J134426.41+441620.0','J131305.81+012755.9','J121826.72-000750.1',
#        'J104252.93+041441.1','J012159.81-010224.3','J084143.50+013149.8','J095540.47+050236.5']
# IDs.sort()
# ##Prepare sample.cat
# bands = ['G', 'R', 'I', 'Z', 'Y']

# mag = []
# text_temp = "# id F314 E314 F315 E315 F316 E316 F317 E317 F318 E318\n"
# for ID in IDs:
#     mag = []
#     idx = np.where(Reines_t1[:,2] == ID)[0][0]
#     Reines_iMag = float(Reines_t1[idx, 5])
#     z = float(Reines_t1[idx, 4])
#     for band in bands:
#         try:
#             fit_run = pickle.load(open(glob.glob('./galight_results/'+ID+'_{0}*pkl'.format(band))[0],'rb')) 
#             obs_mag = fit_run.final_result_galaxy[0]['magnitude']
#         except:
#             obs_mag = -99
#         mag.append(obs_mag)
#     mag_err = [0.1] * len(mag)
#     fnu = [10 ** ((mag[i]-25)/(-2.5)) for i in range(len(mag))]
#     fnu_up = [10 ** ((mag[i]-mag_err[i]-25)/(-2.5)) for i in range(len(mag))]
#     fnu_dw = [10 ** ((mag[i]+mag_err[i]-25)/(-2.5)) for i in range(len(mag))]
#     fnu_err = [(fnu_up[i]-fnu_dw[i])/2 for i in range(len(mag))]
    
#     for i in range(5):
#         if mag[i] <0:
#             fnu[i] = fnu[3]
#             fnu_err[i] = 1000000
    
#     folder_path = 'esti_smass/' + ID + '/'
    # print(ID, z)
    # print(round(fnu[0],3), round(fnu_err[0],3), round(fnu[1],3), round(fnu_err[1],3), round(fnu[2],3), round(fnu_err[2],3),
    #       round(fnu[3],3), round(fnu_err[3],3), round(fnu[4],3), round(fnu_err[4],3)) 
    # os.mkdir(path = folder_path)
#     write_file = open(folder_path+'sample.cat','w') 
#     write_file.write(text_temp)
#     _string = str(int(ID[1:7])) + " {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(round(fnu[0],3), round(fnu_err[0],3), round(fnu[1],3), round(fnu_err[1],3), round(fnu[2],3), round(fnu_err[2],3),
#           round(fnu[3],3), round(fnu_err[3],3), round(fnu[4],3), round(fnu_err[4],3)) 
#     write_file.write(_string)
#     write_file.close()
#     # print(Reines_t1[Reines_t1[:,2]==ID] )  #The stellar mass by Reiness


# for ID in IDs:
#     f = open("temp_esti_stellar_mass/sample.input","r")
#     string = f.read()
#     mag = []
#     idx = np.where(Reines_t1[:,2] == ID)[0][0]
#     Reines_iMag = float(Reines_t1[idx, 5])
#     z = float(Reines_t1[idx, 4])
#     print(ID, z)
#     string = string.replace("1234", str(int(ID[1:7])))
#     string = string.replace("0.0339", str(z))
#     folder_path = 'esti_smass/' + ID + '/'
#     write_file = open(folder_path+'sample.input','w') 
#     write_file.write(string)
#     write_file.close()

# f = open("temp_esti_stellar_mass/run_gsf.py","r")
# string = f.read()
# for ID in IDs:
#     folder_path = 'esti_smass/' + ID + '/'
#     write_file = open(folder_path+'run_gsf.py','w') 
#     write_file.write(string)
#     write_file.close()

# for ID in IDs:
#     print("runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/{0}/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/{0}')".format(ID))
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J012159.81-010224.3/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J012159.81-010224.3')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J084143.50+013149.8/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J084143.50+013149.8')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J095540.47+050236.5/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J095540.47+050236.5')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J104252.93+041441.1/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J104252.93+041441.1')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J121826.72-000750.1/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J121826.72-000750.1')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J131305.81+012755.9/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J131305.81+012755.9')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J131310.12+051942.1/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J131310.12+051942.1')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J134426.41+441620.0/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J134426.41+441620.0')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J140018.41+050242.2/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J140018.41+050242.2')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J141630.82+013708.0/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J141630.82+013708.0')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J141920.64+043623.3/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J141920.64+043623.3')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J143450.63+033842.5/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J143450.63+033842.5')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J145819.55+045451.7/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J145819.55+045451.7')
# runfile('/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J233837.09-002810.4/run_gsf.py', wdir='/Users/Dartoon/Astro/Projects/my_code/projects/2022_HSC_compare_Reines/esti_smass/J233837.09-002810.4')
    