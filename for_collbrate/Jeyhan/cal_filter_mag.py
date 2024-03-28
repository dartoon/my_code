#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:16:20 2024

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob


import sys
sys.path.insert(0,'../../model_z6_data_id0/')

folder = './Generate_template_z0.7/' #
# /Users/Dartoon/Astro/Projects/my_code/projects/2022_JWST_QSOz6/model_z6_data_id1/for_smass/esti_smass/20230328_freeParam_prior2/run_flag6.py


# folder = folder + '_freeParam' #!!!
# folder = folder + '_freeParam_prior1_Av2' #!!!

# folder = folder + '_freeParam_prior1_Av5' #!!! Used in the paper
# folder = folder + '_freeParam_prior1_Av_test' #!!!

# folder = folder + '_freeParam_addF300M_0.6' #
# folder = folder + '_freeParam_withEmissionLine_high-dust' #
# folder = 'esti_smass/202301152_freeParam_withEmissionLine_wide-dustup3' #!!!


steller_file = glob.glob(folder+'/gsf_spec_*.fits')[0]
hdul = pyfits.open(steller_file)
info_muv = hdul[1].header 
print(info_muv['MUV50'], info_muv['MUV16'], info_muv['MUV84'])

steller_file = glob.glob(folder+'/SFH_*.fits')[0]
hdul = pyfits.open(steller_file)
info1 = hdul[0].header 
print('redshift', float(info1['ZMC_50']))
print('smass', info1['Mstel_50'], info1['Mstel_16'], info1['Mstel_84']) 

# smass_l = round(float(info['Mstel_16']),3) 
# smass = round(float(info['Mstel_50']),3) 
# smass_h = round(float(info['Mstel_84']),3) 

#%%

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import seaborn

import sys
hdul = pyfits.open(steller_file)
info = hdul[0].header 
# print(target_id)
print('redshift', float(info['ZMC_50']))
print('smass', float(info['Mstel_50']) )

sfr = round(10**float(info['SFR_50']),3) 
z = info['Z']
sfr_l = round(10**float(info['SFR_16']),3) 
sfr_h = round(10**float(info['SFR_84']),3) 
smass_l = round(float(info['Mstel_16']),3) 
smass = round(float(info['Mstel_50']),3) 
smass_h = round(float(info['Mstel_84']),3) 
age = round(10**float(info['T_MW_50']),3) 
age_l = round(10**float(info['T_MW_16']),3) 
age_h = round(10**float(info['T_MW_84']),3) 
mel = round(float(info['Z_MW_50']),1) 

print('sfr', sfr)
print('sfr_l', sfr_l)
print('sfr_h', sfr_h)
print('age', age)
print('age_l', age_l)
print('age_h', age_h)
print('AV', float(info['AV_50']) )
print('SSFR', round(float(info['SFR_50']) - float(info['Mstel_50']) ,3) )
# print("open",'esti_smass/'+folder+str(idx), '/' )
print('\n')
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
#%%
spec1d_file = glob.glob(folder+'/gsf_spec*.fits')[0]
spec1d = pyfits.open(spec1d_file)  
name_spec = spec1d[1].columns
unit =  spec1d[1].header['TUNIT3']
table_spec = spec1d[1].data
wave = table_spec['wave_model']
# f_16 = table_spec['f_model_16']
f_50 = table_spec['f_model_50']
# f_84 = table_spec['f_model_84']

seaborn.reset_orig()
plt.figure(figsize=(10, 6))
plt.rcParams["font.family"] = "sans-serif"
# array_spec[:,2] =  array_spec[:,2]/ array_spec[:,2].max() * 2.65

plt.plot(wave/10000., f_50, 'black', alpha=0.7)
# plt.fill_between(wave/10000., f_16,f_84, color = 'gray',alpha = 0.5)

hst_filt_id = {'F606W': '4', 'F814W':'6', 'F105W':'202', 'F125W':'203', 'F140W':'204', 'F160W':'205'}

jwst_filt_id = {'F115W': '352', 'F150W': '353', 'F200W': '354', 
           'F277W': '355', 'F356W': '356', 'F444W': '357', 'F410M': '362', 'F480M': '365', 'F300M':'359',
           'F250M': '358'}
all_filt_id = hst_filt_id | jwst_filt_id
filt_id = hst_filt_id | jwst_filt_id

ivd = {v: k for k, v in filt_id.items()}
sample_cat_file = glob.glob(folder+'/sample.cat')[0]

f = open(sample_cat_file,"r")
string = f.read()
lines = string.split('\n')   # Split in to \n
line = lines[0]
filt_id = line.split('F')[1:]
filt_id = [filt_id[i].split(' ')[0] for i in range(len(filt_id))]
fnu_s = lines[1].split(' ')[1:]
fnu_s = [float(fnu_s[i]) for i in range(len(fnu_s))]

flambda_list = []
for i, fid in enumerate(filt_id[::-1]):
    f_fil = np.loadtxt('../../template/gsf_temp/filter/{0}.fil'.format(fid))
    lam = np.median(f_fil[1:,1])
    fnu = fnu_s[::-1][2*i+1]
    fnu_err = fnu_s[::-1][2*i]
    if fnu != 0:
        mag = -2.5*np.log10(fnu) + 25
        flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**float(unit.split('e-')[1][:2])
        lam = lam/10000.
        yerr_l = fnu_err/fnu
        plt.scatter(lam, flambda, c='r', zorder=100,  s=180)
        plt.errorbar(lam, flambda, yerr=yerr_l*flambda,fmt='.',color='red',markersize=1, zorder =100, linewidth=2)
    if fnu == 0:
        mag = -2.5*np.log10(fnu_err) + 25
        flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**float(unit.split('e-')[1][:2])
        lam = lam/10000.
        plt.scatter(lam, flambda, c='r', zorder =100, marker="_", s=280)
        plt.arrow(lam, flambda, 0, -flambda*0.6, length_includes_head=True,
              head_width=0.2, head_length=flambda*0.1, zorder=102, color='red', linewidth=1.2)
        
    f_array = np.vstack((wave, f_50)).T
    filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
    plt.scatter(lam, filt_flam, marker="d", zorder =90,  s=280, facecolors='none', edgecolors='blue', linewidths=2)
    flambda_list.append(flambda)
        


# #plt.plot(sov_jwst_f144w_fil[:,0]/10000., sov_jwst_f144w_fil[:,1], label='NIRCam F444W response curve')
plt.xlim(0.1, 6)
# plt.ylim(0., 0.88)
plt.ylim(0., np.max(flambda_list)*2.0)
xmin, xmax, ymin, ymax = plt.axis()

for i, fid in enumerate(filt_id[::-1]):
    f_fil = np.loadtxt('../../template/gsf_temp/filter/{0}.fil'.format(fid))
    top = f_50.max()
    f_fil[:,2] = f_fil[:,2]/f_fil[:,2].max() * (ymax-ymin) * 0.1
    plt.plot(f_fil[1:,1]/10000., f_fil[1:,2], label='{0}'.format(ivd[fid]))

values = float(info1['Mstel_16']), float(info1['Mstel_50']), float(info1['Mstel_84'])
plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.9, r"log M$_*$="+info['Mstel_50'], fontsize=22)
plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.81, 'age = {0:.1f} Gyr'.format(age, age_l, age_h), fontsize=17)
plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.74, r'A$\mathrm{_{V}}$'+ r'= {0:.1f} '.format(float(info['AV_50'])) , fontsize=17)    
plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.67, r"log Z/Z$_{\rm \odot}$" + r" = {0:.1f}".format(mel), fontsize=17)    
plt.rcParams['legend.edgecolor'] = '0.5'
plt.legend(prop={'size':13}, ncol=2, loc = 1)
plt.tick_params(labelsize=25)
plt.xlabel(r"$\lambda$ ($\mu$m)",fontsize=27)
plt.ylabel(r"f$_\lambda$  (10$^{\rm" + " -{0}}}$".format(unit.split('e-')[1][:2])+" erg s$^{-1}$ cm$^{-2}$$\mathrm{\AA}^{-1}$)",fontsize=27)
# plt.title(target_id,fontsize=27, y=1.02) 
plt.show()

#%% Cal magnitude:
all_filt_id['VIS'] = '5818'
# all_filt_id
str_0 = '# id '
str_1 = '1234 '
filts = ''
list0 = ['F115W', 'F150W', 'F277W', 'F444W', 'VIS']
list1 = ['F814W', 'F115W', 'F150W', 'F277W', 'F444W', 'VIS']
list2 = ['F606W', 'F814W', 'F115W', 'F150W', 'F277W', 'F444W', 'VIS']
for key in all_filt_id.keys():
    # if key in ['F115W', 'F150W', 'F277W', 'F444W']:
    # if key in list0:
    # if key in list1:
    if key in list2:
    # if key in ['F606W','F814W', 'F115W', 'F150W', 'F277W', 'F444W']:
        fid = all_filt_id[key]
        print(key)
        # print(fid)
        f_fil = np.loadtxt('../../template/gsf_temp/filter/{0}.fil'.format(fid))
        filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
        # print(filt_flam)   
        lam = np.median(f_fil[1:,1])
        mag = -2.5* np.log10(filt_flam  / 10**float(unit.split('e-')[1][:2])) - 2.402 - 5.0*np.log10(lam)
        mag_err = 0.3
        if key == 'F814W' and 'F606W' not in list1:
            mag_err = 0.6
        if key == 'F606W':
            mag_err = 0.6
        flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**float(unit.split('e-')[1][:2])
        # print(key, mag)
        fnu = 10 ** ((mag-25)/(-2.5))
        fnu_up = 10 ** ((mag-mag_err-25)/(-2.5))
        fnu_dw = 10 ** ((mag+mag_err-25)/(-2.5))
        fnu_err = (fnu_up-fnu_dw)/2
        # print(fid, fnu, fnu_err)
        str_0 +=  'F{0} E{0} '.format(fid)
        str_1 += '{0:.3f} {1:.3f} '.format(fnu, fnu_err)
        filts += fid + ','
print(str_0[:-1])
print(str_1[:-1])
print(filts[:-1])