#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 13:49:21 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import os, glob
from gsf import gsf
import shutil 

def esti_smass(ID, mags, z, folder = 'esti_smass/'):
    ID = ID
    z = z
    #%%Create a cat file
    folder_path = folder + ID + '/'
    if glob.glob(folder_path) != []:   
        shutil.rmtree(folder_path)
    os.mkdir(path = folder_path)
    text_temp = "# id F314 E314 F315 E315 F316 E316 F317 E317 F318 E318\n"
    mag_err = [0.05] * len(mags)
    fnu = [10 ** ((mags[i]-25)/(-2.5)) for i in range(len(mags))]
    fnu_up = [10 ** ((mags[i]-mag_err[i]-25)/(-2.5)) for i in range(len(mags))]
    fnu_dw = [10 ** ((mags[i]+mag_err[i]-25)/(-2.5)) for i in range(len(mags))]
    fnu_err = [(fnu_up[i]-fnu_dw[i])/2 for i in range(len(mags))]
    for i in range(5):
        if mags[i] <0:
            fnu[i] = 100
            fnu_err[i] = 1000000
    write_file = open(folder_path+'sample.cat','w') 
    write_file.write(text_temp)
    _string = str(int(ID)) + " {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(round(fnu[0],3), round(fnu_err[0],3), round(fnu[1],3), round(fnu_err[1],3), round(fnu[2],3), round(fnu_err[2],3),
          round(fnu[3],3), round(fnu_err[3],3), round(fnu[4],3), round(fnu_err[4],3)) 
    write_file.write(_string)
    write_file.close()
    
    #%%Create a input file
    f = open("../../template/gsf_HSCgrizy_temp/sample.input","r")
    string = f.read()
    string = string.replace("idname", str(int(ID)))
    string = string.replace("zinfo", str(z))
    string = string.replace("folder/", folder_path)
    write_file = open(folder_path+'sample.input','w') 
    write_file.write(string)
    write_file.close()
    gsf.run_gsf_all(folder+'{0}/sample.input'.format(ID), 0, idman=None)
    #%%
    mv_files = glob.glob('*_{0}_*'.format(ID)) + glob.glob('*_{0}.*'.format(ID))
    for mv_file in mv_files:
        shutil.move(mv_file, folder+'{0}/'.format(ID))
    
