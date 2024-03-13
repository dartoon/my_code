#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 20:25:37 2024

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
#Load fit_run  and Rerun

idx = 1  #!!!
import sys
sys.path.insert(0, '../../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

z_str = str(z)

# filters = ['F150W', 'F356W']
filters = ['F356W', 'F150W', 'F115W', 'F250M', 'F444W', 'F200W', 'F300M', 'F480M']
rerun = False
if rerun == True:
    for top_psf_id in range(1):
        for count in range(len(filters)):
            fit_run_list = []
            # idx = idx_info
            filt = filters[count]
            fit_files = glob.glob('../{1}/stage3_all/*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
            if filt == 'F356W' or filt == 'F150W':
                fit_files = glob.glob('../../../2022_JWST_QSOz6/model_z6_data_id{0}/stage3_all/fit_material/fit_run_fixn1__idx1_{1}*pkl'.format(idx, filt))#+\
            fit_files.sort()
            for i in range(len(fit_files)):
                fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
            chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
            sort_Chisq = chisqs.argsort()  
            print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
            fit_run_fix_n1 = fit_run_list[sort_Chisq[top_psf_id]]
            fit_run_fix_n1.plot_final_qso_fit()

            supersampling_factor = 3
            fit_sepc = FittingSpecify(fit_run_fix_n1.fitting_specify_class.data_process_class)
            fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = supersampling_factor,
                                         point_source_supersampling_factor = 1,
                                         apertures_center_focus=True)
            fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
            fit_sepc.kwargs_params['lens_light_model'][4][0]['R_sersic'] = 1.
            fit_sepc.kwargs_params['lens_light_model'][2][0]['n_sersic'] = 4.
            fit_sepc.plot_fitting_sets()
            fit_run = FittingProcess(fit_sepc, savename = target_id)
            fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
            fit_run.plot_final_qso_fit(target_ID =target_id)
            pickle.dump(fit_run , open('./fix_n4_run/fit_run_{0}_idx{1}_psfsp1_fixn4.pkl'.format(filt, idx), 'wb'))

#%%
#Plot The SED:
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
    
folder = '../../../2024_JWST_Cy2/esti_smass/202401261/'
spec1d_file = glob.glob(folder+'/gsf_spec*.fits')[0]
spec1d = pyfits.open(spec1d_file)  
name_spec = spec1d[1].columns
unit =  spec1d[1].header['TUNIT3']
table_spec = spec1d[1].data
wave = table_spec['wave_model']
f_16 = table_spec['f_model_16']
f_50 = table_spec['f_model_50']
f_84 = table_spec['f_model_84']
import seaborn
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
        plt.arrow(lam, flambda, 0, -flambda*0.6, length_includes_head=True,
              head_width=0.2, head_length=flambda*0.1, zorder=102, color='red', linewidth=1.2)
        
    f_array = np.vstack((wave, f_50)).T
    filt_flam = cal_filt_flam(f_array , f_fil[:,1:])    
    plt.scatter(lam, filt_flam, marker="d", zorder =90,  s=280, facecolors='none', edgecolors='blue', linewidths=2)
    flambda_list.append(flambda)
        


# plt.ylim(0., np.max(flambda_list)*2.0)
xmin, xmax, ymin, ymax = plt.axis()

for filt in ['F356W', 'F150W', 'F115W', 'F250M', 'F444W', 'F200W', 'F300M', 'F480M']:
    filename = './fix_n4_run/fit_run_{0}_idx{1}_psfsp1_fixn4.pkl'.format(filt, idx)
    fit_run = pickle.load(open(filename,'rb'))
    mag = fit_run.final_result_galaxy[0]['magnitude']
    fid = jwst_filt_id[filt]
    f_fil = np.loadtxt('../../../../template/gsf_temp/filter/{0}.fil'.format(fid))
    lam = np.median(f_fil[1:,1])
    flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**float(unit.split('e-')[1][:2])
    lam = lam/10000.
    # print(lam, flambda)
    plt.scatter(lam, flambda, c='blue', zorder =100, s=180)
    print(filt, round(mag,3))


for i, fid in enumerate(filt_id[::-1]):
    f_fil = np.loadtxt('../../../../template/gsf_temp/filter/{0}.fil'.format(fid))
    top = f_50.max()
    f_fil[:,2] = f_fil[:,2]/f_fil[:,2].max() * (ymax-ymin) * 0.1
    plt.plot(f_fil[1:,1]/10000., f_fil[1:,2], label='{0}'.format(ivd[fid]))

plt.xlim(0.25, 6)
plt.ylim(0., 1.2)
steller_file = glob.glob(folder+'/SFH_*.fits')[0]
hdul = pyfits.open(steller_file)
info1 = hdul[0].header 

values = float(info1['Mstel_16']), float(info1['Mstel_50']), float(info1['Mstel_84'])
plt.rcParams['legend.edgecolor'] = '0.5'
plt.legend(prop={'size':13}, ncol=2, loc = 1)
plt.tick_params(labelsize=25)
plt.xlabel(r"$\lambda$ ($\mu$m)",fontsize=27)
plt.ylabel(r"f$_\lambda$  (10$^{\rm" + " -{0}}}$".format(unit.split('e-')[1][:2])+" erg s$^{-1}$ cm$^{-2}$$\mathrm{\AA}^{-1}$)",fontsize=27)

