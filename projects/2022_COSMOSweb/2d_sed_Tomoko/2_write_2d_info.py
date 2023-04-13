#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:32:31 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

f = open('fmos_alma_cosmosweb.cat','r')
string = f.read()
lines = string.split('\n')
lines = [lines[i] for i in range(len(lines)) if 'FMOS_J09' in lines[i]]

for i in range(5):
    target_name, RA, Dec, z, best_mass = lines[i].split(' ')
    name = target_name[7:12]
    z = float(z)
    folder = 'esti_smass/{0}/*'.format(name)
    steller_files = glob.glob(folder+'/SFH_*.fits')
    for file in steller_files:
        hdul = pyfits.open(file)
        count = file.split('SFH_')[1].split('.')[0][len(name[4:]):]
        info = hdul[0].header 
        smass = info['Mstel_50']
        sfr = info['SFR_50']
        m_age = info['T_MW_50']
        l_age = info['T_LW_50']
        AV = info['AV_50']
        # shdul = pyfits.open(file.replace('SFH','gsf_spec'))
        # info_stellar = shdul[1].header 
        # BV50 = info_stellar['BV50']
        # f = open( file.split('SFH')[0]+'gsf_spec_header.txt',"r")
        # string = f.read()
        # BV50 = float(string.split('BV50    =')[1].split('BV84')[0])
        filename = name+'_sed_2d_result_bin2.txt'
        if_file = glob.glob(filename)
        if if_file == []:
            write_file =  open(filename,'w')
            write_file.write("count_i, smass, sfr, m_age, l_age, AV, E(B-V) \n")
        else:
            write_file =  open(filename,'r+') 
            write_file.read()
        write_file.write("{0} {1} {2} {3} {4} {5}".format(count, smass, sfr, m_age, l_age, AV))
        write_file.write("\n")
        write_file.close()
