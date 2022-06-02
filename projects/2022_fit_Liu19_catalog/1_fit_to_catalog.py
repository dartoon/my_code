#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 11:05:52 2022

@author: Dartoon
"""


import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import os

bands = 'GRIZY'  #Band that will be download
f = open("catalog.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
import glob
from galight.data_process import DataProcess

write_file = open('table_mag_frameflux.txt','w') 
write_file.write("#ID, (mag using flux with frame) Gmag, Rmag, Imag, Zmag, Ymag  \n")
for _, line in enumerate(lines):
    ID, Ra, Dec = line.split(' ')
    write_file.write(ID+' ')
    for band in bands:
        glob_files = glob.glob('fit_result/{0}-{1}.pkl'.format(ID,band))
        if glob_files != []:
            fit_run  = pickle.load(open(glob_files[0],'rb'))
            # print(fit_run.final_result_galaxy[0])
            write_file.write('{0:.3f} '.format(fit_run.final_result_galaxy[0]['magnitude']))
        else:
            write_file.write('-99 ')
    write_file.write('\n')
write_file.close()

write_file = open('table_mag_sersicflux.txt','w') 
write_file.write("#ID, (mag using flux with frame) Gmag, Rmag, Imag, Zmag, Ymag  \n")
for _, line in enumerate(lines):
    ID, Ra, Dec = line.split(' ')
    write_file.write(ID+' ')
    for band in bands:
        glob_files = glob.glob('fit_result/{0}-{1}.pkl'.format(ID,band))
        if glob_files != []:
            fit_run  = pickle.load(open(glob_files[0],'rb'))
            sersic_mag = -2.5*np.log10(fit_run.final_result_galaxy[0]['flux_sersic_model']) + 27.0
            # print(fit_run.final_result_galaxy[0])
            write_file.write('{0:.3f} '.format(sersic_mag))
        else:
            write_file.write('-99 ')
    write_file.write('\n')
write_file.close()
