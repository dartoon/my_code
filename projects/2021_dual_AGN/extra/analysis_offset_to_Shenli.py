#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 23:33:03 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pandas as pd
import glob

s_sample = pd.read_csv('../Shenli_data/five_band_color.csv', index_col = 0) 

folder = 'NTT_candidates/'
# folder = 'extra_interesting/'
f = open(folder+"cut_out.txt","r")
string = f.read()
cut_out = string.split('\n')   # Split in to \n

files = glob.glob(folder+'*I.fits')
files.sort()

import math
mismatch_overall = []
ID_issue = []
for band in ['I', 'G', 'R']: #, 'Z', 'Y']:
    mismatch = []
    for i in range(len(files)):
        ID = cut_out[i].split(' ')[0]
        #Load Shenli's fitting
        try:
            l = np.where(s_sample['ID'] == ID)[0][0]
        except:
            # print(ID, "exist in source_list.asc but not in five_band_color.csv")
            continue
        s_AGN_0_pos = np.array([s_sample['qso_RA'][l], s_sample['qso_DEC'][l]])
        s_AGN_1_pos = np.array([s_sample['com_RA'][l], s_sample['com_DEC'][l]])
        s_AGN_mag = np.array( [s_sample['qso_{0}'.format(band.lower())][l], s_sample['com_{0}'.format(band.lower())][l]])
        
        file = glob.glob('NTT_candidates/fit_result/{0}/fit_result_{1}-band.txt'.format(ID, band))
        if (math.isnan( s_AGN_mag[0] ) == True or math.isnan( s_AGN_mag[1] ) == True) and file != []:
            print(ID, 'Shenli did not get the band of', band, 'but, it exist')
            # ID_issue.append(ID)
        elif math.isnan( s_AGN_mag[0] ) == False and math.isnan( s_AGN_mag[1] ) == False and file == []:
            print(ID, 'I did not get the band of', band, 'but, Shenli get')
        elif math.isnan( s_AGN_mag[0] ) == False and math.isnan( s_AGN_mag[1] ) == False:
            f = open(file[0],"r")
            Trust_fitting = 2
            string = f.read()
            lines = string.split('\n')   # Split in to \n
            l0 = [ j for j in range(len(lines)) if 'AGN mag:' in lines[j]]
            AGN_mag = lines[l0[Trust_fitting]]
            AGN_mag = AGN_mag.split(' ')[2:4]
            AGN_mag  = np.array([float(AGN_mag[0]), float(AGN_mag[1])])
            l1 = [ j for j in range(len(lines)) if 'AGN0 position:' in lines[j]]
            AGN_0_pos = np.array([float(lines[l1[Trust_fitting-1]].split('RA: ')[1].split(' ')[0]),
                         float(lines[l1[Trust_fitting-1]].split('DEC: ')[1].split(';')[0]) ])
            AGN_1_pos = np.array([float(lines[l1[Trust_fitting-1]].split('RA: ')[2].split(' ')[0]),
                         float(lines[l1[Trust_fitting-1]].split('DEC: ')[2].split(';')[0]) ])
            offset_0 = np.array([np.sqrt(np.sum((s_AGN_0_pos - AGN_0_pos)**2))*3600, np.sqrt(np.sum((s_AGN_0_pos - AGN_1_pos)**2))*3600])
            offset_1 = np.array([np.sqrt(np.sum((s_AGN_1_pos - AGN_0_pos)**2))*3600, np.sqrt(np.sum((s_AGN_1_pos - AGN_1_pos)**2))*3600])
            if np.min(offset_0) > 0.6 or np.min(offset_1) > 0.6:
                test = 0
                print(ID, 's_AGN0 position could not match', 'AGN0_match:', (np.min(offset_0) < 0.6), 
                      'AGN1_match:',(np.min(offset_1) < 0.6), band)
                ID_issue.append(ID)
            else:
                order = np.array([np.where(offset_0 == offset_0.min())[0][0], np.where(offset_1 == offset_1.min())[0][0] ])
                if order[0] == order[1]:
                    print(ID, "There is a position match problem for ID")
                    ID_issue.append(ID)
                else:
                    mag_offset = [s_AGN_mag[0] - AGN_mag[order[0]], s_AGN_mag[1] - AGN_mag[order[1]] ] 
                    mismatch.append(mag_offset)
                    # if np.max(mag_offset) > 1:
                    #     print(ID, band, 'AGN mag mismatch')
                        # ID_issue.append(ID)
    mismatch_overall.append(mismatch)

mismatch_i = np.array(mismatch_overall[0])
ID_issue  = list(dict.fromkeys(ID_issue))
#%%
import shutil
for band in ['I', 'G', 'R']:
    for ID in ID_issue:
        # copy_f = glob.glob('NTT_candidates/fit_result/{0}/fit_{1}-band_fit2_PSPS+Sersic_*.pdf'.format(ID, band))
        # if copy_f != []:
        #     shutil.copy(copy_f[0], '/Users/Dartoon/Downloads/NTT_issue/'+ID+'_{0}-band_fit2.pdf'.format(band))
        shutil.copy('NTT_candidates/fit_result/{0}/fitting2_used_aper.pdf'.format(ID), 
                    '/Users/Dartoon/Downloads/NTT_issue/'+ID+'_fitting2_used_aper.pdf')