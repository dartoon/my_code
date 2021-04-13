#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 18:21:50 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
import copy
import glob

zlines = np.loadtxt('../_pdfs_2close/DR144.4_short.asc', dtype='str')

from ast import literal_eval
def read_string_list(string):
    """
    Translate a string-list-dict to a list-dict. 
    Not allow array inside the dict...
    """
    string = ''.join(string.split('array(['))
    string = ''.join(string.split('])'))    
    string = string.replace('nan', '-99')
    string = string.replace('inf', '-99')
    string_list = string.split('{')[1:]
    string_list = [literal_eval("{"+string_list[i].split('}')[0]+"}") for i in range(len(string_list))]
    return string_list


def read_z(ID):
    line = zlines[zlines[:,0] == ID]
    if len(line) >0:
        z = float(line[0][-1])
    elif len(line) == 0:
        RA_Dec = copy.deepcopy(zlines[:,1:3])
        RA_Dec =  RA_Dec.astype(np.float64)
        pos = SkyCoord('{0} {1}'.format(ID[:2]+':'+ID[2:4]+':'+ID[4:9], ID[9:12]+':'+ID[12:14]+':'+ID[14:]), unit=(u.hourangle, u.deg))
        ID_RA_dec = np.array([pos.ra.degree, pos.dec.degree])
        dis_sq = np.sum((ID_RA_dec - RA_Dec)**2,axis = 1)
        idx = np.where(dis_sq == dis_sq.min())[0][0]
        if dis_sq.min()>(0.35/3600)**2:
            raise ValueError("The searched ID have too much big distance, over 0.35 arcsec: "+ repr(np.sqrt(dis_sq.min()) * 3600 ))
        line  = zlines[idx]
        z = float(line[-1])
    return z

def read_info(ID):
    line = zlines[zlines[:,0] == ID]
    if len(line) >0:
        RA = float(line[0][-3])
        Dec = float(line[0][-2])
        z = float(line[0][-1])
    elif len(line) == 0:
        RA_Dec = copy.deepcopy(zlines[:,1:3])
        RA_Dec =  RA_Dec.astype(np.float64)
        pos = SkyCoord('{0} {1}'.format(ID[:2]+':'+ID[2:4]+':'+ID[4:9], ID[9:12]+':'+ID[12:14]+':'+ID[14:]), unit=(u.hourangle, u.deg))
        ID_RA_dec = np.array([pos.ra.degree, pos.dec.degree])
        dis_sq = np.sum((ID_RA_dec - RA_Dec)**2,axis = 1)
        idx = np.where(dis_sq == dis_sq.min())[0][0]
        line  = zlines[idx]
        line = [x for x in line if x!='']
        if dis_sq.min()>(0.35/3600)**2:
            raise ValueError("The searched ID have too much big distance, over 0.35 arcsec: "+ repr( np.sqrt(dis_sq.min()) * 3600 ))
        RA = float(line[-3])
        Dec = float(line[-2])
        z = float(line[-1])
    return RA, Dec, z

# ID = '001459.72+002319.2'
def cal_oreination(ID):
    files_1 = glob.glob('../proof2close_HSC_images_5band/*/' + ID + '/fit_result/')
    files_2 = glob.glob('../extra/*/fit_result*/' + ID + '/')
    files = files_1 + files_2
    file = glob.glob(files[-1]+'fit_result_{0}-band.txt'.format('I'))
    if file != []:
        f = open(file[0],"r")    
    string = f.read()    
    lines = string.split('\n')   # Split in to \n
    trust = 2    
    l1 = [i for i in range(len(lines)) if 'model_PS_result:' in lines[i]]
    AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])    
    if AGN_dic[0]['flux_within_frame'] < AGN_dic[1]['flux_within_frame']:
        AGN_dic[0], AGN_dic[1] = AGN_dic[1], AGN_dic[0]
    AGN_pos = np.array([[1*AGN_dic[i]['ra_image'], AGN_dic[i]['dec_image']] for i in range(len(AGN_dic))])    
    dif = AGN_pos[1]-AGN_pos[0]
    PA = np.arctan( dif[0]/dif[1] ) * 180 / np.pi
    if dif[1]<0 and dif[0]>0:
        PA = 180 + PA
    if dif[1]<0 and dif[0]<0:
        PA = 180 + PA
    if dif[1]>0 and dif[0]<0:
        PA = 360 + PA
    return PA

# ID = '222057.44+000329.8'
# ID = '020053.43-042721.8'
def cal_offset(ID):
    trust = 2
    offset = None
    file_all_cand_0 = glob.glob('../proofBHBH/allband_data/fit_result/'+ID+'/')
    file_all_cand_1 = glob.glob('../proof2close_HSC_images_5band/z_over1/'+ID+'/')
    file_all_cand =  file_all_cand_0 + file_all_cand_1
    if file_all_cand != []:
        folder = file_all_cand[0]
        if 'proof2close_HSC_images_5band' in folder:
            folder = folder + 'fit_result' + '/'
        file = glob.glob(folder + 'fit_result_{0}-band.txt'.format('I'))    
        if file!=[]:
            f = open(file[0],"r")    
            string = f.read()    
            lines = string.split('\n')   # Split in to \n
            l1 = [i for i in range(len(lines)) if 'model_PS_result:' in lines[i]]
            AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])    
            AGN_pos = np.array([[-1*AGN_dic[i]['ra_image'], AGN_dic[i]['dec_image']] for i in range(len(AGN_dic))])    
            offset = np.sum( (AGN_pos[0] -  AGN_pos[1])**2)**0.5  
    return offset
# print( read_z('162501.98+430931.6') )
# print( read_z('023245.68-000426.1') )
# print(cal_oreination('162501.98+430931.6') )