#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 21:20:09 2017

@author: dxh
"""
import glob

#For change the name from Y7000P, the PSF no error
# names = 'AGN_result_folder/*.pkl'
# folder_list = glob.glob(names)
# folder_list.sort()
# import os
# for i in range(len(folder_list)):
#     ID = folder_list[i].split('/')[1][:3]
#     name = folder_list[i]
#     name1 = folder_list[i].split('/')[0] + '/' + 'idx{0}_ID{1}_PSFnoerr_PSFinter.pkl'.format(i, ID)
#     if glob.glob(name) != []:
#         os.rename(name, name1)

# #For copy the previous PSF var 0.01
# from shutil import copyfile
# names = 'simulations_700_subg30/sim_lens_ID_subg30_*/result_PSFerr001_PSFinter_subg3.pkl'
# folder_list = glob.glob(names)
# folder_list.sort()
# for i in range(len(folder_list)):
#     name = folder_list[i]
#     ID = name.split('/')[1][-3:]
#     name1 = 'AGN_result_folder/' + 'idx{0}_ID{1}_PSFerr001_PSFinter.pkl'.format(i, ID)
#     if glob.glob(name) != []:
#         copyfile(name, name1)
