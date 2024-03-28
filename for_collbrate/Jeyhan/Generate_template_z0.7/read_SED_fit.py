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
sys.path.insert(0,'../../model_z6_data_id0/')

folder = './' #
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

#%%
import pickle, glob
file = glob.glob(folder+'/chain_*_phys.cpkl')
value = pickle.load(open(file[0],'rb'))
import pandas as pd 

seaborn.set(font_scale=1.5)

prop = list(value['chain'].keys())
for i in range(len(prop)):
    if i == 1:
        value['chain'][prop[i]] = 10**value['chain'][prop[i]]

df = pd.DataFrame(data=value['chain'])
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
for key in df.keys():
    values.append([np.percentile(df[key], 16), np.median(df[key]), np.percentile(df[key], 84)])
g = seaborn.pairplot(df,corner=True, diag_kind="kde", plot_kws=dict(marker="+", s=20, linewidth=1))
# g.map_lower(seaborn.kdeplot, levels=4, color=".2")
    
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
            # title = title.replace('0.00', '0.01')
            ax.set_title(title, fontsize=20)
            # ax.set_title(df.keys()[i/5]+r'$\displaystyle\substack{1\2}$'.format(values[j][1], values[j][2]-values[j][1], values[j][1]-values[j][0]))
            # ax.set_title(r'\frac{-e^{i\pi}}{2^n}$!')
            # ax.set_title(r'\frac{-e^{i\pi}}{2^n}$!', fontsize=16, color='r')
            # ax.set_title(r'\TeX\ is Number $\displaystyle\sum_{n=1}^\infty'
            #              r'\frac{-e^{i\pi}}{2^n}$!', fontsize=16, color='r')
            # ax.set_title(df.keys()[i/5]+r'={0}$\pm$')
        # if j %5 ==0:
            # ax.set_xlim(10.2, 11.9)
        if Mstr == '':
            Mstr = title
# plt.savefig('../../model_z6_data_id0/figures/{0}_SED_MCMC.pdf'.format(target_id[:5]), bbox_inches = "tight")

#%%

# import matplotlib as mat
# mat.rcParams['font.family'] = 'STIXGeneral'

# folder = '202209' #0.2 mag no HST.
# folder = '20220901_' #Not HST
# folder = '20220901' #HST upper limit
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
spec1d_file = glob.glob(folder+'/gsf_spec*.fits')[0]
spec1d = pyfits.open(spec1d_file)  
name_spec = spec1d[1].columns
unit =  spec1d[1].header['TUNIT3']
table_spec = spec1d[1].data
wave = table_spec['wave_model']
f_16 = table_spec['f_model_16']
f_50 = table_spec['f_model_50']
f_84 = table_spec['f_model_84']

seaborn.reset_orig()
plt.figure(figsize=(10, 6))
plt.rcParams["font.family"] = "sans-serif"
# array_spec[:,2] =  array_spec[:,2]/ array_spec[:,2].max() * 2.65

plt.plot(wave/10000., f_50, 'black', alpha=0.7)
plt.fill_between(wave/10000., f_16,f_84, color = 'gray',alpha = 0.5)

hst_filt_id = {'F606W': '4', 'F814W':'6', 'F105W':'202', 'F125W':'203', 'F140W':'204', 'F160W':'205'}

jwst_filt_id = {'F115W': '352', 'F150W': '353', 'F200W': '354', 
           'F277W': '355', 'F356W': '356', 'F444W': '357', 'F410M': '362', 'F480M': '365', 'F300M':'359',
           'F250M': '358'}
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
    f_fil = np.loadtxt('../../../template/gsf_temp/filter/{0}.fil'.format(fid))
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
    # f_array = np.vstack((wave, f_16)).T
    # filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
    # plt.scatter(lam, filt_flam, marker="d", zorder =90,  s=280, facecolors='none', edgecolors='green', linewidths=2)
    # f_array = np.vstack((wave, f_84)).T
    # filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
    # plt.scatter(lam, filt_flam, marker="d", zorder =90,  s=280, facecolors='none', edgecolors='green', linewidths=2)
    
    flambda_list.append(flambda)
        


# #plt.plot(sov_jwst_f144w_fil[:,0]/10000., sov_jwst_f144w_fil[:,1], label='NIRCam F444W response curve')
plt.xlim(0.1, 6)
plt.ylim(0., 0.88)
# plt.ylim(0., np.max(flambda_list)*2.0)
xmin, xmax, ymin, ymax = plt.axis()

for i, fid in enumerate(filt_id[::-1]):
    f_fil = np.loadtxt('../../../template/gsf_temp/filter/{0}.fil'.format(fid))
    top = f_50.max()
    f_fil[:,2] = f_fil[:,2]/f_fil[:,2].max() * (ymax-ymin) * 0.1
    plt.plot(f_fil[1:,1]/10000., f_fil[1:,2], label='{0}'.format(ivd[fid]))


# plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.93, 'z: {0}'.format(z), fontsize=17)
# plt.text( (xmax-xmin)*005, (ymax-ymin)*0.87, r"M$_*$: [{1:.2f}$\leftarrow${0:.2f}$\rightarrow${2:.2f}]".format(smass, smass_l, smass_h), fontsize=17)
# plt.text( (xmax-xmin)*0.05, (ymax-ymin)*0.80, "M$_{uv}$: "+ r"[{1:.2f}$\leftarrow${0:.2f}$\rightarrow${2:.2f}]".format(info_muv['MUV50'], info_muv['MUV84'], info_muv['MUV16']), fontsize=17)
# print('M$_*$:  [{1:.2f}, {0:.2f} ,{2:.2f}]'.format(smass, smass_l, smass_h))
# print('Muv:  [{1:.2f}, {0:.2f} ,{2:.2f}]'.format(info_muv['MUV50'], info_muv['MUV84'], info_muv['MUV16']))
# if '{0:.1f}'.format(age_l) == '0.0':
#     plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.61, r'age: [$\leftarrow${2:.1f}]  Gyr'.format(age, age_l, age_h), fontsize=17)
#     print(r'age:  [$\leftarrow${2:.1f}]  Gyr'.format(age, age_l, age_h))
# else:
#     plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.61, 'age fixed as {0:.1f} Gyr'.format(age, age_l, age_h), fontsize=17)
#     print('age:  [{1:.1f}$-${2:.1f}]  Gyr'.format(age, age_l, age_h))
    # print('age:  [{1:.1f}$-${2:.1f}]  Gyr'.format(age, age_l, age_h))
values = float(info1['Mstel_16']), float(info1['Mstel_50']), float(info1['Mstel_84'])
plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.9, r"log M$_*$="+Mstr, fontsize=22)
# plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.95, r"log M$_*$ = {0:.2f}".format(smass), fontsize=17)
# plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.88, "M$\mathrm{_{uv}}$ = "+ r"{0:.2f}".format(info_muv['MUV50']), fontsize=17)
plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.81, 'age = {0:.1f} Gyr'.format(age, age_l, age_h), fontsize=17)
plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.74, r'A$\mathrm{_{V}}$'+ r'= {0:.1f} '.format(float(info['AV_50'])) , fontsize=17)    
plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.67, r"log Z/Z$_{\rm \odot}$" + r" = {0:.1f}".format(mel), fontsize=17)    
# plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.81, 'age = {0:.1f} Gyr'.format(values[1][1]), fontsize=17)
# plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.74, r'A$\mathrm{_{V}}$'+ r'= {0:.1f}'.format(float(info['AV_50'])) , fontsize=17)    
# plt.text( (xmax-xmin)*0.07, (ymax-ymin)*0.67, r"log Z/Z$_{\rm \odot}$" + r" = {0:.1f}".format(mel), fontsize=17)    
plt.rcParams['legend.edgecolor'] = '0.5'
plt.legend(prop={'size':13}, ncol=2, loc = 1)
plt.tick_params(labelsize=25)
plt.xlabel(r"$\lambda$ ($\mu$m)",fontsize=27)
plt.ylabel(r"f$_\lambda$  (10$^{\rm" + " -{0}}}$".format(unit.split('e-')[1][:2])+" erg s$^{-1}$ cm$^{-2}$$\mathrm{\AA}^{-1}$)",fontsize=27)
# plt.title(target_id,fontsize=27, y=1.02) 
# plt.savefig('../../model_z6_data_id0/figures/{0}_SED_map.pdf'.format(target_id[:5]), bbox_inches = "tight")
