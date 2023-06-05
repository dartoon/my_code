#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 23:48:59 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob


import sys
sys.path.insert(0,'../')

folder = '20230328' #
fitidx = 0
from target_info import target_info
info = target_info[str(fitidx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

idx = 1  #our target 102 fix age.

folder = 'esti_smass/'+folder #+str(idx)
folder = folder + '_freeParam_withEmissionLine_Av0-1'

steller_file = glob.glob(folder+'/gsf_spec_*.fits')[0]
hdul = pyfits.open(steller_file)
info_muv = hdul[1].header 

steller_file = glob.glob(folder+'/SFH_*.fits')[0]
hdul = pyfits.open(steller_file)
info1 = hdul[0].header 


#%%

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob,seaborn

#%%
import pickle, glob
import pandas as pd 
import copy

seaborn.set(font_scale=1.5)
value_list = []
for i in range(5):
    file = glob.glob('esti_smass/20230328_freeParam_withEmissionLine_Av{0}-{1}/chain_*_phys.cpkl'.format(i,i+1))
    value = pickle.load(open(file[0],'rb'))
    prop = list(value['chain'].keys())
    for i in range(len(prop)):
        if i == 1:
            value['chain'][prop[i]] = 10**value['chain'][prop[i]]
    value_list.append(value)
    
    
value = copy.deepcopy(value_list[0]['chain'])
for key in value.keys():
    for i in range(1,5):
        value[key] = np.concatenate([value[key],value_list[i]['chain'][key]])

bool_array = np.random.lognormal(0, 0.8, size=value['AV'].shape) > value['AV']
for key in value.keys():
    value[key] = value[key][bool_array]

df = pd.DataFrame(data=value)
plt.rcParams["font.family"] = "sans-serif"
for key in df.keys():
    if key == 'logM_stel':
        newkey = r"log M$_{*}$/M$_{\rm \odot}$"
        df[newkey] = df.pop(key)
    if key == 'logZ_MW':
        newkey = r"log Z/Z$_{\rm \odot}$"
        df[newkey] = df.pop(key)
    if key == 'AV':
        newkey = r'A$\mathrm{_{V}}$'
        df[newkey] = df.pop(key)
    if key == 'logT_MW':
        newkey = 'age (Gyr)'
        df[newkey] = df.pop(key)
values = []
for i,key in enumerate(df.keys()):
    if i != 0:
        values.append([np.percentile(df[key], 16), np.median(df[key]), np.percentile(df[key], 84)])
    if i == 0 :
        values.append([10.53-0.37, 10.53, 10.53+0.51])
g = seaborn.pairplot(df,corner=True, diag_kind="kde", plot_kws=dict(marker="+", s=20, linewidth=1))

Mstr = ''   
for i in range(len(g.axes.ravel())):
    ax = g.axes.ravel()[i]
    if ax != None:
        ax.spines['top'].set_visible(True) 
        ax.spines['right'].set_visible(True)
        ax.spines['left'].set_visible(True)
        j = i % len(values)
        ax.axvline(x=values[j][1], ls='--', linewidth=1.6, c='coral')
        ax.axvline(x=values[j][0], ls='--', linewidth=1.6, c='coral')
        ax.axvline(x=values[j][2], ls='--', linewidth=1.6, c='coral')
        ax.set_xlabel(ax.get_xlabel(), fontsize=25)
        ax.set_ylabel(ax.get_ylabel(), fontsize=25)
        ax.tick_params(labelsize=18) 
        if i %5 == 0:
            title_fmt=".2f"
            fmt = "{{0:{0}}}".format(title_fmt).format
            title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
            title = title.format(fmt(values[j][1]), fmt(values[j][1]-values[j][0]), fmt(values[j][2]-values[j][1]))
            # if j == 0:
            #     title = '${10.53}_{-0.37}^{+0.51}$'
            ax.set_title(title, fontsize=20)
            # title = title.replace('0.00', '0.01')
        if j %5 ==0:
            ax.set_xlim(9.4, 11.8)
        if Mstr == '':
            Mstr = title
            # Mstr = '${10.53}_{-0.37}^{+0.51}$'
plt.savefig('../../model_z6_data_id0/figures/{0}_SED_MCMC.pdf'.format(target_id[:5]), bbox_inches = "tight")


# import matplotlib as mat
# mat.rcParams['font.family'] = 'STIXGeneral'

#%%

import sys
sys.path.insert(0,'../')

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
f_16_l, f_50_l, f_84_l = [], [], []
for i in range(2):
    spec1d_file = glob.glob('esti_smass/20230328_freeParam_withEmissionLine_Av{0}-{1}'.format(i,i+1)+'/gsf_spec*.fits')[0]
    spec1d = pyfits.open(spec1d_file)  
    name_spec = spec1d[1].columns
    unit =  spec1d[1].header['TUNIT3']
    table_spec = spec1d[1].data
    wave = table_spec['wave_model']
    f_16_l.append(table_spec['f_model_16'])
    f_50_l.append(table_spec['f_model_50'])
    f_84_l.append(table_spec['f_model_84'])
    #%%
f_16 = np.min(f_16_l, axis = 0)
f_84 = np.max(f_84_l, axis = 0)
# test_i = 1
# f_16 = f_16_l[test_i]
# f_84 = f_84_l[test_i]
f_50 = f_50_l[0]
    

seaborn.reset_orig()
plt.figure(figsize=(10, 6))
plt.rcParams["font.family"] = "sans-serif"
# array_spec[:,2] =  array_spec[:,2]/ array_spec[:,2].max() * 2.65

plt.plot(wave/10000., f_50, 'black', alpha=0.7)
plt.fill_between(wave/10000., f_16,f_84, color = 'gray',alpha = 0.5)

hst_filt_id = {'F606W': '4', 'F814W':'6', 'F105W':'202', 'F125W':'203', 'F140W':'204', 'F160W':'205'}

jwst_filt_id = {'F115W': '352', 'F150W': '353', 'F200W': '354', 
           'F277W': '355', 'F356W': '356', 'F444W': '357', 'F410M': '362'}
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
        plt.errorbar(lam, flambda, yerr=yerr_l*flambda,fmt='.',color='red',markersize=1, zorder =100, linewidth=2)
    if fnu == 0:
        mag = -2.5*np.log10(fnu_err) + 25
        flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**float(unit.split('e-')[1][:2])
        lam = lam/10000.
        plt.scatter(lam, flambda, c='r', zorder =100, marker="_", s=280)
        plt.arrow(lam, flambda, 0, -flambda*0.2, length_includes_head=True,
              head_width=0.2, head_length=flambda*0.1, zorder=102, color='red', linewidth=1.2)
        
        
    f_array = np.vstack((wave, f_50)).T
    filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
    plt.scatter(lam, filt_flam, marker="d", zorder =90,  s=280, facecolors='none', edgecolors='blue', linewidths=2)
    # f_array = np.vstack((wave, f_16)).T
    # filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
    # plt.scatter(lam, filt_flam, marker="d", zorder =90,  s=280, facecolors='none', edgecolors='green', linewidths=2)
    # f_array = np.vstack((wave, f_84)).T
    # filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
    # plt.scatter(lam, filt_flam, marker="d", zorder =90,  s=280, facecolors='none', edgecolors='green', linewidths=2)

    # flambda_list.append(flambda)
        
        
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
plt.ylim(0., 0.31)
# plt.ylim(0., np.max(flambda_list)*2.0)
xmin, xmax, ymin, ymax = plt.axis()

for i, fid in enumerate(filt_id):
    f_fil = np.loadtxt('../../../../template/gsf_temp/filter/{0}.fil'.format(fid))
    top = f_50.max()
    f_fil[:,2] = f_fil[:,2]/f_fil[:,2].max() * (ymax-ymin) * 0.1
    plt.plot(f_fil[1:,1]/10000., f_fil[1:,2], label='{0}'.format(ivd[fid]))
values = float(info1['Mstel_16']), float(info1['Mstel_50']), float(info1['Mstel_84'])
plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.9, r"log M$_*$ = "+Mstr, fontsize=22)
plt.rcParams['legend.edgecolor'] = '0.5'
plt.legend(prop={'size':20}, ncol=2, loc = 1)
plt.tick_params(labelsize=25)
plt.xlabel(r"$\lambda$ ($\mu$m)",fontsize=27)
plt.ylabel(r"f$_\lambda$  (10$^{\rm" + " -{0}}}$".format(unit.split('e-')[1][:2])+" erg s$^{-1}$ cm$^{-2}$$\mathrm{\AA}^{-1}$)",fontsize=27)
plt.title(target_id,fontsize=27, y=1.02) 
plt.savefig('../figures/{0}_SED_map.pdf'.format(target_id[:5]), bbox_inches = "tight")
#plt.yticks([])
