#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 17:02:15 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

#%%Check the z>1 ones that have 2close and have proof_BHBH
folders = ['233718.07+002550.6', '230402.77-003855.4', '222929.45+010438.4', '221101.45+001449.0', 
          '220642.82+003016.2', '162501.98+430931.6', '153008.91+425634.8', '150216.66+025719.8', 
          '135944.20+012809.8', '134257.16-013912.9', '133222.62+034739.9', '132441.58-015401.8', 
          '124618.51-001750.2', #'113803.73+031457.8', 
          '105458.01+043310.6', '104644.31+000329.7', 
          '101138.50+012344.6', '090654.53+021315.2', '022906.04-051428.9', '022404.85+014941.9',
          '020318.87-062321.3', '020231.14-042246.8', '020053.43-042721.8', '015141.69-000646.8', 
          '014235.36+012334.0', '013834.18-000509.3', '013736.57+005742.3', '013526.15+004135.8',
          '011227.87-003151.6', '001401.62-002600.5']
for ID in folders:
    path = '../_John_fitted/' + ID + '*/proof-BHBH*'
    file = glob.glob(path)
    print(ID, file)
    
#%%Check for the z>1 ones that have proof_BHBH:
files = glob.glob('../_John_fitted/*')

f = open("../_pdfs_2close/DR144.4_short.asc","r")
string = f.read()
lines = string.split('\n')   # Split in to \n

def read_z(ID):
    line = [lines[i] for i in range(len(lines)) if ID in lines[i]]
    if line != []:
        z = float(line[0].split(' ')[-1])
    else:
        z = -99
    return z

proof_BHBHs = []
for file in files:
    ID = (file.split('/')[-1]).split('_HSC')[0]
    if read_z(ID)>1:
        proof_BHBH = glob.glob('../_John_fitted/' + ID + '*/proof-BHBH*')
        if proof_BHBH!= []:
            proof_BHBHs.append(proof_BHBH[0])
#%%
# import shutil
# for file in proof_BHBHs:
#     copyfile = '/Users/Dartoon/Downloads/proof_BHBH/' + file.split('/')[-2]+'-'+file.split('/')[-1] + '.pdf'
#     shutil.copy(file, copyfile)        
            