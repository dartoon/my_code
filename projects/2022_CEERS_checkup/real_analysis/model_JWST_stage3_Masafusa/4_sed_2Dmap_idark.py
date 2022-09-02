#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 09:48:51 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
from functions_for_result import esti_smass, load_prop, load_info
import glob
import os # Delete xfile.txt
import matplotlib
matplotlib.use('Agg')
import sys
count = int(sys.argv[1]) - 1 # 1 - 2500
flag = int(sys.argv[2])

idx = 1
target_id, z = load_info(idx)
sed_2d_info = pickle.load(open('sed_2d_info_bin2.pkl','rb'))
mag_dict = sed_2d_info[count][2]
esti_smass(ID = '20220903'+str(int(count)), mags_dict = mag_dict, z = z, flag = flag, if_run_gsf=True)
folder = 'esti_smass/20220903'+str(count)
spec_file = glob.glob(folder+'/gsf_spec_*.fits')[0]
hdul_spec = pyfits.open(spec_file)
info_spec = hdul_spec[1].header

write_file = open(folder + '/gsf_spec_header.txt','w') 
write_file.write(str(info_spec))
write_file.close()
rm_file = glob.glob(folder+'/gsf_spec_*.fits') + glob.glob(folder + '/SPEC*corner*png') + glob.glob(folder+'/*asdf') 
for file in rm_file:
    os.remove(file)
# steller_file = glob.glob(folder+'/SFH_*.fits')[0]
# hdul = pyfits.open(steller_file)
# info = hdul[0].header 
# smass = info['Mstel_50']
# sfr = info['SFR_50']
# m_age = info['T_MW_50']
# l_age = info['T_LW_50']
# AV = info['AV_50']
# filename = 'sed_2d_result_bin2.txt'
# if_file = glob.glob(filename)
# if if_file == []:
#     write_file =  open(filename,'w')
#     write_file.write("count_i, smass, sfr, m_age, l_age, AV \n")
# else:
#     write_file =  open(filename,'r+') 
#     write_file.read()
# write_file.write("{0} {1} {2} {3} {4} {5}".format(count, smass, sfr, m_age, l_age, AV))
# write_file.write("\n")
# write_file.close()