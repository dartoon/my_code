#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 10:50:03 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pandas as pd
import glob


folder = 'NTT_candidates/'
# folder = 'extra_interesting/'
f = open(folder+"cut_out.txt","r")
string = f.read()
cut_out = string.split('\n')   # Split in to \n

files = glob.glob(folder+'*I.fits')
files.sort()

from ast import literal_eval
def read_string_list(string):
    """
    Translate a string-list-dict to a list-dict. 
    Not allow array inside the dict...
    """
    string = ''.join(string.split('array(['))
    string = ''.join(string.split('])'))    
    string = string.replace('nan', '-99')
    string_list = string.split('{')[1:]
    string_list = [literal_eval("{"+string_list[i].split('}')[0]+"}") for i in range(len(string_list))]
    return string_list
    
    

import math
mismatch_overall = []
ID_issue = []
for i in range(len(files)):
    for band in ['I', 'G', 'R']: #, 'Z', 'Y']:
        ID = cut_out[i].split(' ')[0]
        #Load Shenli's fitting
        file = glob.glob('NTT_candidates/fit_result/{0}/fit_result_{1}-band.txt'.format(ID, band))
        if file != []:
            f = open(file[0],"r")
            Trust_fitting = 2
            string = f.read()
            lines = string.split('\n')   # Split in to \n
            # l0 = [ j for j in range(len(lines)) if 'AGN mag:' in lines[j]]
            AGN_dic = read_string_list(string = lines[-3].split('model_PS_result: ')[1])
            AGN_mag = [AGN_dic[i]['magnitude'] for i in range(2)]
            AGN_0_pos = np.array([AGN_dic[0]['ra_image'], AGN_dic[0]['dec_image']])
            AGN_1_pos = np.array([AGN_dic[1]['ra_image'], AGN_dic[1]['dec_image']])
            galaxy_dic = read_string_list(string = lines[-4].split('model_Sersic_result: ')[1])
            galaxy_mag = [galaxy_dic[i]['magnitude'] for i in range(len(galaxy_dic))]
            galaxy_pos = np.array([[galaxy_dic[i]['center_x'], galaxy_dic[i]['center_y']] for i in range(len(galaxy_dic))])
            
            # AGN_mag = lines[l0[Trust_fitting]]
            # AGN_mag = AGN_mag.split(' ')[2:4]
            # AGN_mag  = np.array([float(AGN_mag[0]), float(AGN_mag[1])])
            # l1 = [ j for j in range(len(lines)) if 'AGN0 position:' in lines[j]]
            # AGN_0_pos = np.array([float(lines[l1[Trust_fitting-1]].split('RA: ')[1].split(' ')[0]),
            #              float(lines[l1[Trust_fitting-1]].split('DEC: ')[1].split(';')[0]) ])
            # AGN_1_pos = np.array([float(lines[l1[Trust_fitting-1]].split('RA: ')[2].split(' ')[0]),
            #              float(lines[l1[Trust_fitting-1]].split('DEC: ')[2].split(';')[0]) ])
            distance_0 = np.sqrt(np.sum((AGN_0_pos - galaxy_pos)**2, axis = 1))
            galaxy_ID_0 = np.where(distance_0 == distance_0.min())[0][0]
            distance_1 = np.sqrt(np.sum((AGN_1_pos - galaxy_pos)**2, axis = 1))
            galaxy_ID_1 = np.where(distance_1 == distance_1.min())[0][0]            
            
            # if distance < 1:
            #     print(ID, band, 'too close')
            
