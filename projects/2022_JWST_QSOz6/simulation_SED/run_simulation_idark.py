#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:59:54 2023

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import shutil,os,glob
from scipy.interpolate import interp1d
def cal_filt_flam(array_spec, fil):
    '''
    Calculate a filter f_lambda
    Parameter
    --------
    Return
    --------
        Filter flux
    '''
    fil[:,1] /= fil[:,1].max()
    f_filt = interp1d(fil[:,0], fil[:,1], kind='cubic')
    int_spec = array_spec[(array_spec[:,0]>fil[0,0]) * (array_spec[:,0]<fil[-1,0])]
    int_flux = 0
    int_filt = 0
    for i in range(len(int_spec)-1):
        int_flux += int_spec[i,1] * f_filt(int_spec[i,0]) * (int_spec[i+1,0]-int_spec[i,0])
        int_filt += f_filt(int_spec[i,0]) * (int_spec[i+1,0]-int_spec[i,0])
    filt_flam = int_flux/ int_filt   
    return filt_flam

#%%Produce mock paramters

import sys
count = int(sys.argv[1]) - 1 # 1 - 2809
filter_set = int(sys.argv[3])
# count = 1
# filter_set = 0

# flag = int(sys.argv[2])
flag = 1
if filter_set == 0:
    filts = ['F150W', 'F356W']
if filter_set == 1:
    filts = ['F277W', 'F356W']
if filter_set == 2:
    filts = ['F356W', 'F444W']
if filter_set == 3:
    filts = ['F150W', 'F277W', 'F356W']
if filter_set == 4:
    filts = ['F150W', 'F277W', 'F356W', 'F444W']
if filter_set == 5:
    filts = ['F150W', 'F277W']
if filter_set == 6:
    filts = ['F150W', 'F200W']
    
#%%
seed = count
age = np.random.uniform(0.3,0.7)#0.3
Av = np.random.uniform(0.3,1.0)
# metallicity = np.random.uniform(-0.7, -0.3)

folder = 'second_run/seed{0}_sim'.format(seed)
if filter_set == 0:
    if glob.glob(folder) !=[]:
        shutil.rmtree(folder)
    os.mkdir(folder)
    # shutil.copy('gene_temp/sample.cat', folder)
    f = open("gene_temp/sample.cat","r")
    string = f.read()
    string = string.replace("seedid", str(seed) )
    write_file = open(folder+'/generate_sample.cat','w') 
    write_file.write(string)
    write_file.close()
    
    #Create a input file
    f = open("gene_temp/sample.input","r")
    string = f.read()
    string = string.replace("age_temp", str(age))
    string = string.replace("folder", folder)
    # string = string.replace("Z_temp", str(metallicity))
    string = string.replace("Av_temp", str(Av) )
    string = string.replace("seedid", str(seed) )
    write_file = open(folder+'/generate_sample.input','w') 
    write_file.write(string)
    write_file.close()
    
    from gsf import gsf
    gsf.run_gsf_all(folder + '/generate_sample.input', 0, idman=None)
    move_file = glob.glob('*_20230301{0}_*'.format(seed)) + glob.glob('*_20230301{0}.*'.format(seed))
    for file in move_file:
        shutil.move(file, folder)
#%%calculate mags for SED fit
import glob
spec1d_file = glob.glob(folder+'/gsf_spec_*.fits')[0]
spec1d = pyfits.open(spec1d_file)  
name_spec = spec1d[1].columns
unit =  spec1d[1].header['TUNIT3']
table_spec = spec1d[1].data
steller_file = glob.glob(folder+'/SFH_*.fits')[0]
hdul = pyfits.open(steller_file)
info1 = hdul[0].header 
smass_True = info1['Mstel_50']
z = info1['Z']
wave = table_spec['wave_model']
f_50 = table_spec['f_model_50']
jwst_filt_id = {'F115W': '352', 'F150W': '353', 'F200W': '354', 
            'F277W': '355', 'F356W': '356', 'F444W': '357', 'F410M': '362'}
filters = ''
mag_result =  {}#{'F150W': 26.4, 'F356W': 24.80} #QSO host
for filt in filts:
    fid = jwst_filt_id[filt]
    f_fil = np.loadtxt('../../../template/gsf_temp/filter/{0}.fil'.format(fid))
    lam = np.median(f_fil[1:,1])
    f_array = np.vstack((wave, f_50)).T
    filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
    flambda = filt_flam 
    mag= -2.5*np.log10(flambda / (10**float(unit.split('e-')[1][:2]))) - (2.402+5.0*np.log10(lam)) 
    print(filt, mag)
    mag_result[filt] = mag + np.random.normal(0,0.2) #!!! Adding uncertainty
    filters = filters + filt + '_'
filters = filters[:-1]    
#%%SED fit
#Create a input file
SED_folder = folder.split('/')[0] + '/' + filters + '/'
if glob.glob(SED_folder) ==[]:
    os.mkdir(SED_folder)

from SED_function import esti_smass
sed_folder = SED_folder + folder.split('/')[1].replace('sim', 'result')

esti_smass(ID = '20230302'+str(seed),folder = folder.replace('sim', 'result'), 
           mags_dict = mag_result, z = z, flag = flag, 
            if_run_gsf=True, band_as_upper = [],
            mag_err=[0.2]*len(mag_result), just_run = False)
#%%Collecting results
steller_file = glob.glob(folder.replace('sim', 'result')+'/SFH_*.fits')[0]
hdul = pyfits.open(steller_file)
info1 = hdul[0].header 
smass_infer = info1['Mstel_50']