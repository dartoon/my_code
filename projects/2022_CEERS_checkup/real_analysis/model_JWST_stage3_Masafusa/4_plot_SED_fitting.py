#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 00:03:18 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from functions_for_result import esti_smass, load_prop, load_info, name_list

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

folder = '20220901'


idx = 2

_, z = load_info(idx)
target_id = name_list[idx]

print(idx)
steller_file = glob.glob('esti_smass/'+folder+str(idx)+'/SFH_*.fits')[0]
hdul = pyfits.open(steller_file)
info = hdul[0].header 
print(target_id)
print('redshift', float(info['ZMC_50']))
print('smass', float(info['Mstel_50']) )

sfr = round(10**float(info['SFR_50']),3) 
sfr_l = round(10**float(info['SFR_16']),3) 
sfr_h = round(10**float(info['SFR_84']),3) 
age = round(10**float(info['T_MW_50']),3) 
age_l = round(10**float(info['T_MW_16']),3) 
age_h = round(10**float(info['T_MW_84']),3) 

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
#%%
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

spec1d_file = glob.glob('esti_smass/'+folder+str(idx)+'/gsf_spec*.fits')[0]
spec1d = pyfits.open(spec1d_file)  
name_spec = spec1d[1].columns
unit =  spec1d[1].header['TUNIT3']
table_spec = spec1d[1].data
wave = table_spec['wave_model']
f_16 = table_spec['f_model_16']
f_50 = table_spec['f_model_50']
f_84 = table_spec['f_model_84']

plt.figure(figsize=(10, 6))

# array_spec[:,2] =  array_spec[:,2]/ array_spec[:,2].max() * 2.65

plt.plot(wave/10000., f_16, 'gray', alpha=0.4)
plt.plot(wave/10000., f_50, 'black', alpha=0.7)
plt.plot(wave/10000., f_84, 'gray', alpha=0.4)

hst_filt_id = {'F606W': '4', 'F814W':'6', 'F105W':'202', 'F125W':'203', 'F140W':'204', 'F160W':'205'}

jwst_filt_id = {'F115W': '352', 'F150W': '353', 'F200W': '354', 
           'F277W': '355', 'F356W': '356', 'F444W': '357', 'F410M': '362'}
filt_id = hst_filt_id | jwst_filt_id

ivd = {v: k for k, v in filt_id.items()}
sample_cat_file = glob.glob('esti_smass/'+folder+str(idx)+'/sample.cat')[0]

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
    f_fil = np.loadtxt('../../../../template/gsf_temp/filter/{0}.fil'.format(fid))
    lam = np.median(f_fil[1:,1])
    fnu = fnu_s[::-1][2*i+1]
    fnu_err = fnu_s[::-1][2*i]
    if fnu != 0:
        mag = -2.5*np.log10(fnu) + 25
        flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**float(unit.split('e-')[1][:2])
        lam = lam/10000.
        yerr_l = fnu_err/fnu
        plt.scatter(lam, flambda, c='r', zorder=100,  s=180)
        plt.errorbar(lam, flambda, yerr=yerr_l*flambda,fmt='.',color='gray',markersize=1, zorder =90, linewidth=2)
    if fnu == 0:
        mag = -2.5*np.log10(fnu_err) + 25
        flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**float(unit.split('e-')[1][:2])
        lam = lam/10000.
        plt.scatter(lam, flambda, c='r', zorder =100, marker="_", s=280)
        plt.arrow(lam, flambda, 0, -flambda*0.6, length_includes_head=True,
              head_width=0.2, head_length=flambda*0.1, zorder=102, color='red', linewidth=1.2)
        
    f_array = np.vstack((wave, f_50)).T
    filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
    plt.scatter(lam, filt_flam, marker="d", zorder =98,  s=280, facecolors='none', edgecolors='blue', linewidths=2)
    
    flambda_list.append(flambda)
        
        
# fnu = np.array([316.14340, 329.82858, 17.06660, 44.70776])
# lam = np.array([13971.05, 15418.99, 5373.10, 8084.25])
# yerr_l = np.array([155.291/562.887, 155.291/562.887, 155.291/562.887, 155.291/562.887])
# mag = -2.5*np.log10(fnu) + 25
# flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**17
# lam = lam/10000.
# plt.scatter(lam, flambda, c='r', zorder =100)
# plt.errorbar(lam, flambda, yerr=yerr_l*flambda,fmt='.',color='gray',markersize=1, zorder =90)


# #plt.plot(sov_jwst_f144w_fil[:,0]/10000., sov_jwst_f144w_fil[:,1], label='NIRCam F444W response curve')
plt.xlim(0.25, 6)
plt.ylim(0., np.max(flambda_list)*1.6)
xmin, xmax, ymin, ymax = plt.axis()

for i, fid in enumerate(filt_id[::-1]):
    f_fil = np.loadtxt('../../../../template/gsf_temp/filter/{0}.fil'.format(fid))
    top = f_50.max()
    f_fil[:,2] = f_fil[:,2]/f_fil[:,2].max() * (ymax-ymin) * 0.1
    plt.plot(f_fil[1:,1]/10000., f_fil[1:,2], label='{0} response'.format(ivd[fid]))

# plt.text( (xmax-xmin)*0.7, (ymax-ymin)*0.53, 'stellar population with:', fontsize=17)
plt.text( (xmax-xmin)*0.6, (ymax-ymin)*0.53, 'z={0}'.format(z), fontsize=17)
plt.text( (xmax-xmin)*0.6, (ymax-ymin)*0.43, 'age=[{1:.3f}, {0:.3f}, {2:.3f}]  Gyr'.format(age, age_l, age_h), fontsize=17)
plt.text( (xmax-xmin)*0.6, (ymax-ymin)*0.33, "SFR=[{1:.3f} {0:.3f} {2:.3f}]".format(sfr, sfr_l, sfr_h) + r" M$_{\rm \odot}$/yr ", fontsize=17)
# plt.text( (xmax-xmin)*0.6, (ymax-ymin)*0.25, r"M$_{\rm \odot}$/yr ", fontsize=17)
# plt.text( (xmax-xmin)*0.7, (ymax-ymin)*0.25, r"$24.13\substack{+1.17\\\\-0.55}$", fontsize=17)
# plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.37, 'Z$_*$/Z$_\odot$={0}'.format(mel), fontsize=17)
plt.legend(prop={'size':15}, ncol=2, loc = 1)
plt.tick_params(labelsize=20)
plt.xlabel("um",fontsize=27)
plt.ylabel(r"f$_\lambda$ 10$^{\rm" + " -{0}}}$ erg/s/cm$^2$/$\AA$".format(unit.split('e-')[1][:2]),fontsize=27)
plt.title(target_id,fontsize=27)
plt.savefig('outcomes/idx{0}_sed.pdf'.format(idx))
# #plt.yticks([])


