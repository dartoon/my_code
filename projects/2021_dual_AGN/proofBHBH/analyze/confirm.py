#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 11:45:21 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import shutil
folders =  glob.glob('../model_Iband/*')
folders.sort()

AGN_mag_list = []
AGN_pos_list = []
ID_list = []
offset_list = []
for folder in folders:
    # folder = folders[1]    
    file = glob.glob(folder + '/fit_result_*txt')[0]
    f = open(file,"r")
    string = f.read()
    lines = string.split('\n')   # Split in to \n
    l0 = [ j for j in range(len(lines)) if 'AGN mag:' in lines[j]]
    AGN_mag = lines[l0[2]]
    AGN_mag = AGN_mag.split(' ')[2:4]
    AGN_mag = [float(AGN_mag[0]), float(AGN_mag[1])]
    # AGN_mag_list.append(AGN_mag)
    l1 = [ j for j in range(len(lines)) if 'PS PS' in lines[j]]
    offset = lines[l1[1]].split(' ')[-1]
    l2 = [ j for j in range(len(lines)) if 'model_PS_result:' in lines[j]]
    pos_ = lines[l2[2]]
    x0 = float(pos_.split("'ra_image': array([")[1].split('])')[0])
    x1 = float(pos_.split("'ra_image': array([")[2].split('])')[0])
    y0 = float(pos_.split("'dec_image': array([")[1].split('])')[0])
    y1 = float(pos_.split("'dec_image': array([")[2].split('])')[0])
    AGN_pos = [[-1*x0/0.168, y0/0.168], [-1*x1/0.168, y1/0.168]]
    AGN_pos = np.array(AGN_pos)
    
    l3 = [ j for j in range(len(lines)) if 'Comp x:' in lines[j]]
    comp_x = lines[l3[0]][8:].split(', ')
    l4 = [ j for j in range(len(lines)) if 'Comp y:' in lines[j]]
    comp_y = lines[l4[0]][8:].split(', ')
    comp_pos = []
    
    l5 = [ j for j in range(len(lines)) if 'Comp mag:' in lines[j]]
    comp_mag = lines[l5[0]][10:].split(', ')
    comp_mag = [float(comp_mag[i]) for i in range(len(comp_mag))]
    
    for i in range(len(comp_x)):
        comp_pos.append([-1*float(comp_x[i])/0.168, float(comp_y[i])/0.168])
    comp_pos = np.array(comp_pos)    
    
    confirm = False
    deltaMag = []
    for i in range(len(AGN_pos)):
        dis = np.sqrt(np.sum((AGN_pos[i] - comp_pos)**2, axis = 1))
        if dis.min() < 4:
            idx = np.where(dis == dis.min())[0][0]
            deltaMag.append(AGN_mag[i] - comp_mag[idx])
        else:
            deltaMag.append(-99)
            
    # offset_list.append(float(offset))
    if np.max(AGN_mag) < 22.5 and np.max(deltaMag) < 1:
        print(file)
        ID = file.split('model_Iband/')[1].split('/')[0]
        ID_list.append(ID)
        file1 = glob.glob(file.split('fit_result_I-band.txt')[0] + 'fit_I-band_fit2*')[0]
        # if float(offset) < 1:
        #     shutil.copy(file1, '/Users/Dartoon/Downloads/proofBHBH_interesting/'+ID+file1.split('-band_')[-1])
        # offset_list.append(float(offset))
    else:
        if float(offset) < 1:
            shutil.copy(file1, '/Users/Dartoon/Downloads/proofBHBH_interesting/_'+ID+file1.split('-band_')[-1])                
        
        
        # file2 = glob.glob(ID_list + 'fitting2_used_aper.pdf')[0]
        # shutil.copy(file2, '/Users/Dartoon/Downloads/proofBHBH_interesting/'+ID+file2.split('/')[-1])