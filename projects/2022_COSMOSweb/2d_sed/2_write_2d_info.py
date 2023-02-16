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

for name in ['cid_473', 'cid_1210', 'cid_1245']:
    folder = 'esti_smass/{0}/*'.format(name)
    steller_files = glob.glob(folder+'/SFH_*.fits')
    for file in steller_files:
        hdul = pyfits.open(file)
        count = file.split('/')[2][4:]
        info = hdul[0].header 
        smass = info['Mstel_50']
        sfr = info['SFR_50']
        m_age = info['T_MW_50']
        l_age = info['T_LW_50']
        AV = info['AV_50']
        # shdul = pyfits.open(file.replace('SFH','gsf_spec'))
        # info_stellar = shdul[1].header 
        f = open( file.split('SFH')[0]+'gsf_spec_header.txt',"r")
        string = f.read()
        BV50 = float(string.split('BV50    =')[1].split('BV84')[0])
        filename = name+'_sed_2d_result_bin2.txt'
        if_file = glob.glob(filename)
        if if_file == []:
            write_file =  open(filename,'w')
            write_file.write("count_i, smass, sfr, m_age, l_age, AV, E(B-V) \n")
        else:
            write_file =  open(filename,'r+') 
            write_file.read()
        write_file.write("{0} {1} {2} {3} {4} {5} {6:.4f}".format(count, smass, sfr, m_age, l_age, AV, BV50))
        write_file.write("\n")
        write_file.close()
