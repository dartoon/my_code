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
import shutil

files = glob.glob('../z_over1/*')

for file in files:
    fit_result = file + '/fit_result/' + 'fit_result_I-band.txt'
    fit_result = glob.glob(fit_result)
    if fit_result != []:
        ID = fit_result[0].split('z_over1/')[1].split('/')[0]
        f = open(fit_result[0],"r")
        string = f.read()
        lines = string.split('\n')   # Split in to \n    
        l1 = [ j for j in range(len(lines)) if 'PS PS' in lines[j]]
        if len(l1) > 1:
            offset = lines[l1[1]].split(' ')[-1]
            print(offset)
            fit_image = file + '/fit_result/' + 'fit_I-band_fit2_PSPS+Sersic*'
            fit_image  = glob.glob(fit_image)[0]
            print(fit_image)
            copy_to = '/Users/Dartoon/Downloads/proof_2close_interesting/'+ID+'_I-band_fit2_PSPS+Sersic.pdf'
            print(copy_to )
            if float(offset ) < 1:
                shutil.copy(fit_image, copy_to)