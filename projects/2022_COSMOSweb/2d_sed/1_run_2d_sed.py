#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 17:53:43 2023

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import glob
import os # Delete xfile.txt
import matplotlib
matplotlib.use('Agg')
import sys
count = int(sys.argv[1]) - 1 # 1 - 2809
flag = int(sys.argv[2])
target_ID = int(sys.argv[3])
# count = 4001 # 1 - 2809
# flag = 1
# target_ID = 1

cata_list = pickle.load(open('../material/cata_list.pkl','rb'))
if target_ID == 1:
    check_name= 'cid_473'  #29
if target_ID == 2:
    check_name= 'cid_1210' #8
if target_ID == 3:
    check_name= 'cid_1245' #10
check_id = [i for i in range(len(cata_list)) if cata_list[i][-1] == check_name]
# print(cata_list[check_id[0]])
z = cata_list[check_id[0]][-2]
sed_2d_info = pickle.load(open('2d_filts_mag_bin2_{0}.pkl'.format(check_name),'rb'))

mag_dict = sed_2d_info[count][2]

#%%
from functions_for_result import esti_smass
folder = 'esti_smass/'+check_name+'/'+check_name[4:]+str(count)
if glob.glob(folder+'/SFH_*.fits') == [] and len(mag_dict)>2:
    esti_smass(ID = check_name[4:]+str(int(count)), folder = 'esti_smass/'+check_name+'/',
                mags_dict = mag_dict, z = z, flag = flag, if_run_gsf=True)
    spec_file = glob.glob(folder+'/gsf_spec_*.fits')[0]
    hdul_spec = pyfits.open(spec_file)
    info_spec = hdul_spec[1].header
    write_file = open(folder + '/gsf_spec_header.txt','w') 
    write_file.write(str(info_spec))
    write_file.close()
    rm_file = glob.glob(folder+'/gsf_spec_*.fits') + glob.glob(folder + '/SPEC*corner*png') + glob.glob(folder+'/*asdf') 
    for file in rm_file:
        os.remove(file)