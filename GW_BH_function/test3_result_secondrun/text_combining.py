#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 15:59:21 2019

@author: Dartoon

Combining all the text together
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob
filenames = 'test3_mode1_take2_level20_p*.txt'
files = glob.glob(filenames)
files.sort()
texts = []
for file_i in files:
    with open(file_i) as f:
        content = f.readlines()
    texts.append(content)
text = []

filename = 'test3_mode1_take2_level20.txt'
result = open(filename,'w') 
for i in range(len(texts)):
    for j in range(len(texts[i])):
#        text.append(texts[i][j])
        result.write(texts[i][j])    
result.close()