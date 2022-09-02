#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 11:47:15 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
write_file =  open('sed_2d_result_bin4.txt','w')
write_file.write("count_i, smass, sfr, m_age, l_age, AV \n")
steller_files = glob.glob('esti_smass/20220830*/SFH_*.fits')

steller_files.sort()
for steller_file in steller_files:
    hdul = pyfits.open(steller_file)
    count = steller_file.split('20220830')[1].split('/')[0]
    info = hdul[0].header 
    smass = info['Mstel_50']
    sfr = info['SFR_50']
    m_age = info['T_MW_50']
    l_age = info['T_LW_50']
    AV = info['AV_50']
    write_file.write("{0} {1} {2} {3} {4} {5}".format(count, smass, sfr, m_age, l_age, AV))
    write_file.write("\n")
write_file.close()