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
count = int(sys.argv[1]) - 1 # 1 - 8100
flag = int(sys.argv[2]) #0-1
target_ID = int(sys.argv[3]) #0-5
# count = 4001 # 1 - 2809
# flag = 1
# target_ID = 0

f = open('fmos_alma_cosmosweb.cat','r')
string = f.read()
lines = string.split('\n')
lines = [lines[i] for i in range(len(lines)) if 'FMOS_J09' in lines[i]]

cata_list = pickle.load(open('../material/cata_list.pkl','rb'))

target_name, RA, Dec, z, best_mass = lines[target_ID].split(' ')
z = float(z)

t_name = target_name[5:12]
sed_2d_info = pickle.load(open('2d_filts_mag_bin2_{0}.pkl'.format(t_name),'rb'))
mag_dict = sed_2d_info[count][2]

#%%
name = t_name[2:]
from functions_for_result import esti_smass
folder = 'esti_smass/'+name+'/'+name+str(count)
if glob.glob('esti_smass/'+name) == []:
    os.mkdir(path = 'esti_smass/'+name)

if glob.glob(folder+'/SFH_*.fits') == [] and len(mag_dict)>2:
    esti_smass(ID = name+str(int(count)), folder = 'esti_smass/'+name+'/',
                mags_dict = mag_dict, z = z, flag = flag, band_as_upper = [],if_run_gsf=True)
    spec_file = glob.glob(folder+'/gsf_spec_*.fits')[0]
    hdul_spec = pyfits.open(spec_file)
    info_spec = hdul_spec[1].header
    write_file = open(folder + '/gsf_spec_header.txt','w') 
    write_file.write(str(info_spec))
    write_file.close()
    rm_file = glob.glob(folder+'/gsf_spec_*.fits') + glob.glob(folder + '/*png') + glob.glob(folder+'/*asdf')  + glob.glob(folder+'/*cpkl') 
    for file in rm_file:
        os.remove(file)
        
rm_file = glob.glob('./templates/*{0}.asdf'.format(name[1:]+str(int(count))))
for file in rm_file:
    os.remove(file)