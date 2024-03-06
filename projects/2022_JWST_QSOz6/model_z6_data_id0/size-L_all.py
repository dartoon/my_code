#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 20:48:49 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib as matt
import glob
matt.rcParams['font.family'] = 'STIXGeneral'


#shibuya2015, lg(r0/kpc) =-0.38+ -0.56,  r0 = 0.419+-1.12  beta = 0.25+-0.05

plt.figure(figsize=(11, 7))
filt_wave_dic = {'F115W': 11623.89, 'F150W': 15104.23, 'F200W': 20028.15, 
           'F277W': 27844.64, 'F356W': 35934.49, 'F444W': 44393.50, 'F410M': 40886.55}
wave_rest = 1600

mass_list = []
Reff_list = []

# #Target 0
# Muv_0 = [-19.47, -18.75, -18.05]
# R_filt_0 = np.array([1.87, 1.06 ])
# q = 0.630
# R_filt_0 = R_filt_0*np.sqrt(q)
# # wave = F356W
# smass = 10.63
# z = 6.34
# key = 'F356W'
# filt_wave = filt_wave_dic[key]
# wave = filt_wave / (1+z)
# z_p = int(key[1:4])/(wave_rest/100)  -1
# Reff_filt = R_filt_0[0]
# vdep = -0.35 + 0.12 * z - 0.25 * np.log10( 10 ** smass / 10**10)
# R_0 = [Reff_filt * ((1+z)/(1+z_p)) ** vdep]
# R_0.append(R_0[0]* R_filt_0[1]/R_filt_0[0] )
# plt.scatter(Muv_0[0], np.log10(R_0[0]),c='blue',s=480,marker="o",
#             zorder=100, edgecolors='black',alpha=1, label = 'J2255+0251')
# plt.errorbar(Muv_0[0], np.log10(R_0[0]), 
#              yerr = [ [np.log10(R_0[0]) - np.log10(R_0[0]-R_0[1] )], 
#                      [ np.log10(R_0[0]+R_0[1] ) - np.log10(R_0[0])] ],
#                xerr = [ [Muv_0[1] - Muv_0[0]], [Muv_0[2]-Muv_0[1]] ],
#               c='blue', elinewidth=3, zorder =100)

# #Target 1
# Muv_0 = [-20.97, -20.86, -20.62]
# key = 'F150W'
# # key = 'F356W'
# if key == 'F356W':
#     R_filt_0 = np.array([0.75, 0.11])
#     q = 0.360
# if key == 'F150W':
#     R_filt_0 = np.array([0.59, 0.14])
#     q = 0.269
# R_filt_0 = R_filt_0*np.sqrt(q)
# smass = 11.14
# z = 6.4
# filt_wave = filt_wave_dic[key]
# wave = filt_wave / (1+z)
# z_p = int(key[1:4])/(wave_rest/100)  -1
# Reff_filt = R_filt_0[0]
# vdep = -0.35 + 0.12 * z - 0.25 * np.log10( 10 ** smass / 10**10)
# R_0 = [Reff_filt * ((1+z)/(1+z_p)) ** vdep]
# R_0.append(R_0[0]* R_filt_0[1]/R_filt_0[0] )
# plt.scatter(Muv_0[0], np.log10(R_0[0]),c='orange',s=480,marker="o",
#             zorder=100, edgecolors='black',alpha=1, label = 'J2236+0032')
# plt.errorbar(Muv_0[0], np.log10(R_0[0]), 
#               yerr = [ [np.log10(R_0[0]) - np.log10(R_0[0]-R_0[1] )], 
#                       [ np.log10(R_0[0]+R_0[1] ) - np.log10(R_0[0])] ],
#                 # xerr = [ [Muv_0[1] - Muv_0[0]], [Muv_0[2]-Muv_0[1]] ],
#                 xerr = [ [0.15], [0.2] ],
#                 c='orange', zorder=100,
#               elinewidth=3)
import pickle 
from target_info import target_info
for idx in range(10):
    # Muv_0 = [-20.97, -20.86, -20.62]
    # smass = 11.14
    # z = 6.4
    info = target_info[str(idx)]
    target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
    folder = '../esti_stellar_mass/gsf_id_{0}/'.format(idx)
    steller_file = glob.glob(folder+'/gsf_spec_*.fits')[0]
    hdul = pyfits.open(steller_file)
    info_muv = hdul[1].header 
    
    steller_file = glob.glob(folder+'/SFH_*.fits')[0]
    hdul = pyfits.open(steller_file)
    info1 = hdul[0].header 
    
    smass = float(info1['Mstel_50'])
    Muv_0 = [info_muv['MUV84'], info_muv['MUV50'], info_muv['MUV16']]
    # Muv_0 = [-20.97, -20.86, -20.62]
    
    #%% Load photo
    # idx = idx_info
    top_psf_id = 0
    run_folder = '../model_z6_data_id{0}/stage3_all/'.format(idx) #!!!
    for filt in ['F150W', 'F356W']:
        fit_run_list = []
        key = filt
        PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
        if idx != 1:
            fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        else:
            fit_files = glob.glob(run_folder+'*fit_material/fit_run*_fixn1_*idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        # print('idx', idx, target_id, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]
        deltaPix = fit_run.fitting_specify_class.deltaPix
        fit_run.savename = 'figures/' + fit_run.savename+'_'+filt
        count_n = 5
        fit_files[sort_Chisq[top_psf_id]] 
        Chisq_best = chisqs[sort_Chisq[top_psf_id]]
        Chisq_last= chisqs[sort_Chisq[count_n-1]]
        inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
        weight = np.zeros(len(chisqs))
        for i in sort_Chisq[:count_n]:
            weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
        prop_name = 'flux_within_frame'
        all_values = [100*(fit_run_list[i].final_result_galaxy[0][prop_name]/ (fit_run_list[i].final_result_galaxy[0][prop_name] + fit_run_list[i].final_result_ps[0][prop_name]))
                      for i in range(len(fit_run_list))]
        weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        host_ratio = weighted_value
        rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        print('host flux ratio ' , ": {0:.1f}%+-{1:.1f}%".format(weighted_value, rms_value))    
        if host_ratio > 5 or filt == 'F356W':
        # if weighted_value > 5 or filt == 'F150W':
        # if weighted_value > 5:
            prop_name = 'R_sersic'
            # all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
            all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
            weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
            rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
            # print('host Reff ":' , "{0:.2f}+-{1:.2f}".format(weighted_value, rms_value))
            R_filt_0 = np.array([weighted_value, rms_value])
            
            prop_name = 'q'
            all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
            weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
            rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
            # print('host', prop_name, ": {0:.2f}+-{1:.2f}".format(weighted_value, rms_value))
            q = weighted_value
            break
        
    print(key)
    import astropy.units as u
    from astropy.cosmology import LambdaCDM, FlatLambdaCDM
    cosmo1 = LambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, 0.7)
    arc_per_kpc = cosmo1.arcsec_per_kpc_proper(z).value
    R_filt_0 = R_filt_0/arc_per_kpc
    
    R_filt_0 = R_filt_0*np.sqrt(q)
    filt_wave = filt_wave_dic[key]
    wave = filt_wave / (1+z)
    z_p = int(key[1:4])/(wave_rest/100)  -1
    Reff_filt = R_filt_0[0]
    vdep = -0.35 + 0.12 * z - 0.25 * np.log10( 10 ** smass / 10**10)
    if host_ratio>5:
        R_0 = [Reff_filt * ((1+z)/(1+z_p)) ** vdep]
        R_0.append(R_0[0]* R_filt_0[1]/R_filt_0[0] )
        if idx ==0:
            label = "SHELLQS"
        else:
            label = None
        plt.scatter(Muv_0[0], np.log10(R_0[0]),c='orange',s=480,marker="o",
                    zorder=100, edgecolors='black',alpha=1, label = label)
        plt.errorbar(Muv_0[0], np.log10(R_0[0]), 
                      yerr = [ [np.log10(R_0[0]) - np.log10(R_0[0]-R_0[1] )], 
                              [ np.log10(R_0[0]+R_0[1] ) - np.log10(R_0[0])] ],
                        # xerr = [ [Muv_0[1] - Muv_0[0]], [Muv_0[2]-Muv_0[1]] ],
                        xerr = [ [0.15], [0.2] ],
                        c='orange', zorder=100,
                       elinewidth=3)    
        mass_list.append(smass)
        Reff_list.append(R_0)


#%%
M=np.linspace(-24,-17,10)
r0_s = [-0.2, 0.25]
beta_s = [0.27, 0.01]
r0_liter = [r0_s]
beta_liter = [beta_s]
color_liter = ['black']
label_liter = [' Shibuya+15']
k=0
r0sim=r0_liter[k][0] + r0_liter[k][1] * np.random.randn(100)
bsim=beta_liter[k][0] + beta_liter[k][1] * np.random.randn(100)
X=[]
for i in range(100):
    X.append(-bsim[i]*0.4*(M+21)+r0sim[i])
X=np.array(X)
mu=X.mean(axis=0)
sigma=X.std(axis=0)
plt.plot(M, mu, lw=5, color=color_liter[k], label =label_liter[k],zorder=1)
plt.fill_between(M, (mu+sigma), (mu-sigma), facecolor=color_liter[k], alpha=0.15)
#%%Plot Shibuya:
f = open("../other_files/shibuya","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
lines = [line for line in lines if 'z6' in line or 'z7' in line]

import astropy.units as u
from astropy.cosmology import LambdaCDM, FlatLambdaCDM
cosmo1 = LambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, 0.7)

arc_per_kpc = {6:cosmo1.arcsec_per_kpc_proper(6).value, 7: cosmo1.arcsec_per_kpc_proper(7).value}
distmod = {6:cosmo1.distmod(6).value, 7: cosmo1.distmod(7).value}
# 5*np.log10(cosmo1.luminosity_distance(z=6).value * 10**6) -5

shibuya_scatter = []
for line in lines:
    z = float(line[1])
    info = line.split(' ')
    info = [float(info[i]) for i in range(1,len(info)) if '' != info[i]]
    muv = info[0]
    emuv = info[1]
    Muv = muv - distmod[z] + 2.5 * np.log10(1+z) 
    reff = info[2] * np.sqrt(info[4]) / arc_per_kpc[z] # info[4] is q.
    ereff = info[3] * np.sqrt(info[4])  / arc_per_kpc[z] # info[4] is q.
    shibuya_scatter.append( [Muv, emuv, reff, ereff])
shibuya_scatter = np.array(shibuya_scatter)
plt.scatter(shibuya_scatter[:,0], np.log10(shibuya_scatter[:,2]),c='gray',s=180,marker=".",
            zorder=0, edgecolors='black',alpha=0.6, label = 'Shibuya')
#%%
f = open("../other_files/kawamata","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
lines = [line for line in lines if '+' in line ]
kawamata_scatter = []
for line in lines:
    info = line.split('\t')
    kawamata_scatter.append([float(info[4][:6]), float(info[5][:4])])
kawamata_scatter = np.array(kawamata_scatter)
plt.scatter(kawamata_scatter[:,0], np.log10(kawamata_scatter[:,1]),c='orange',s=180,marker=".",
            zorder=0, edgecolors='black',alpha=0.6, label = 'Kawamata+18')

#%%
sample, idx67, idx_up = np.load('../other_files/HFFsample_yang2022_forDXH.npy',allow_pickle=True)
color_list = ['black']*6
cluster_list=['A2744cl', 'M0416cl', 'M1149cl', 'M0717cl', 'RXC2248', 'A370']
# plt.figure(figsize=(14,9))    
for i in range(6):
    upplinear=(np.log10(sample['Rkpcmcmc'][idx67][sample['cluster'][idx67]==cluster_list[i]]
                   +sample['Rkpc84th'][idx67][sample['cluster'][idx67]==cluster_list[i]])
              -np.log10(sample['Rkpcmcmc'][idx67][sample['cluster'][idx67]==cluster_list[i]])
          )

    lowlinear=(np.log10(sample['Rkpcmcmc'][idx67][sample['cluster'][idx67]==cluster_list[i]])-
                  np.log10(sample['Rkpcmcmc'][idx67][sample['cluster'][idx67]==cluster_list[i]]
                  -sample['Rkpc16th'][idx67][sample['cluster'][idx67]==cluster_list[i]])
          )
    if i ==0:
        label = 'Yang+22'
    else:
        label =None
    plt.scatter((sample['UV_Mag'][idx67][sample['cluster'][idx67]==cluster_list[i]]),
                np.log10(sample['Rkpcmcmc'][idx67][sample['cluster'][idx67]==cluster_list[i]]),
                # xerr=sample['f160w_magerr'][idx67][sample['cluster'][idx67]==cluster_list[i]] * 0,
                # yerr=np.array([ lowlinear, upplinear])*0,
                marker = '.',
                c= 'red',s=180,zorder=1, edgecolors='k', alpha=0.6, label = label)
    # plt.quiver(sample['UV_Mag'][idx_up][sample['cluster'][idx_up]==cluster_list[i]],
    #           np.log10(sample['Rkpcmcmc'][idx_up][sample['cluster'][idx_up]==cluster_list[i]]),
    #           sample['UV_Mag'][idx_up][sample['cluster'][idx_up]==cluster_list[i]]*0,
    #           sample['UV_Mag'][idx_up][sample['cluster'][idx_up]==cluster_list[i]]*0-0.5,
    #           color=color_list[i],
    #            width=0.002,zorder=1, alpha=0.5
    #           )

#%%
plt.xlabel(r"M$_{\rm UV}$",fontsize =30)#,weight='bold')
plt.ylabel(r"$R_{\mathrm{eff, circ.,}}$$_{1600~\AA}$, (log kpc)",fontsize =30)#,weight='bold')
plt.tick_params(labelsize=20)
plt.ylim([-1.5, 1.5])
plt.xlim([-23, -17])
plt.legend(scatterpoints=1,numpoints=1,loc=1,prop={'size':20},ncol=2,handletextpad=0)
# plt.savefig('figures/z6-7_size_luminosity.pdf')
plt.show()


#%%
plt.figure(figsize=(11, 7))
for i in range(len(mass_list)):
    mass =  mass_list[i]
    R_0 = Reff_list[i]
    plt.scatter(mass, np.log10(R_0[0]),c='orange',s=480,marker="o",
                zorder=100, edgecolors='black',alpha=1, label = label)
    plt.errorbar(mass, np.log10(R_0[0]), 
                  yerr = [ [np.log10(R_0[0]) - np.log10(R_0[0]-R_0[1] )], 
                          [ np.log10(R_0[0]+R_0[1] ) - np.log10(R_0[0])] ],
                    # xerr = [ [Muv_0[1] - Muv_0[0]], [Muv_0[2]-Muv_0[1]] ],
                    xerr = [ [0.15], [0.2] ],
                    c='orange', zorder=100,
                   elinewidth=3)    
plt.xlabel(r"M$_{\rm star}$",fontsize =30)#,weight='bold')
plt.ylabel(r"$R_{\mathrm{eff, circ.,}}$$_{1600~\AA}$, (log kpc)",fontsize =30)#,weight='bold')
plt.tick_params(labelsize=20)
plt.ylim([-1.5, 1.5])
plt.xlim([7.6, 11.4])
plt.show()
