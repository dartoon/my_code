#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 15:21:04 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

folders = ['233718.07+002550.6', '230402.77-003855.4', '222929.45+010438.4', '221101.45+001449.0', 
          '220642.82+003016.2', '162501.98+430931.6', '153008.91+425634.8', '150216.66+025719.8', 
          '135944.20+012809.8', '134257.16-013912.9', '133222.62+034739.9', '132441.58-015401.8', 
          '124618.51-001750.2', '113803.73+031457.8', '105458.01+043310.6', '104644.31+000329.7', 
          '101138.50+012344.6', '090654.53+021315.2', '022906.04-051428.9', '022404.85+014941.9',
          '020318.87-062321.3', '020231.14-042246.8', '020053.43-042721.8', '015141.69-000646.8', 
          '014235.36+012334.0', '013834.18-000509.3', '013736.57+005742.3', '013526.15+004135.8',
          '011227.87-003151.6', '001401.62-002600.5']

# move_folder = '../interesting_sources' 
# import shutil
# for i in range(len(folders)):
#     filename = '../z_over1/' + folders[i] + '/fit_result/images_5_band.pdf'
#     shutil.copy(filename, '../interesting_sources/'+ folders[i] +'_images_5_band.pdf')
    
    
folders.sort()

f = open("../../pdfs_2close/DR144.4_short.asc","r")
string = f.read()
lines = string.split('\n')   # Split in to \n

def read_z(ID):
    line = [lines[i] for i in range(len(lines)) if ID in lines[i]]
    if line != []:
        z = float(line[0].split(' ')[-1])
    else: 
        z = -99
    return z    

for ID in folders:
    print(ID, read_z(ID))
    
# ['001401.62-002600.5',
#  '011227.87-003151.6',
#  '013526.15+004135.8',
#  '013736.57+005742.3',
#  '013834.18-000509.3',
#  '014235.36+012334.0',
#  '015141.69-000646.8',
#  '020053.43-042721.8',
#  '020231.14-042246.8',
#  '020318.87-062321.3',
#  '022404.85+014941.9',
#  '022906.04-051428.9',
#  '090654.53+021315.2',
#  '101138.50+012344.6',
#  '104644.31+000329.7',
#  '105458.01+043310.6',
#  '113803.73+031457.8',
#  '124618.51-001750.2',
#  '132441.58-015401.8',
#  '133222.62+034739.9',
#  '134257.16-013912.9',
#  '135944.20+012809.8',
#  '150216.66+025719.8',
#  '153008.91+425634.8',
#  '162501.98+430931.6',
#  '220642.82+003016.2',
#  '221101.45+001449.0',
#  '222929.45+010438.4',
#  '230402.77-003855.4',
#  '233718.07+002550.6'   #color different 
#  ]