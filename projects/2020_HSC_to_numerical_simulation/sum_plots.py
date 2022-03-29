#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:35:10 2021

@author: Dartoon
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sum_prep import load_HSC_comp, load_HST_comp, imf_dict
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

import glob
import pickle
from matplotlib.ticker import AutoMinorLocator
import seaborn as sns


m_ml, b_ml = (0.981139684856507, -2.545890295477823)
#%%
eps = 1.5
fig = plt.figure(figsize=(13*eps, 20*eps))
gs = fig.add_gridspec(6, 4, hspace=0, wspace=0)
axs = gs.subplots()
i, j = 0, 0
obsname = ['HSC', 'HST']
colors_sim = ['deepskyblue', 'steelblue', 'c', 'deeppink', 'hotpink', 'm']

Imag_list = {0.3: 19.5, 0.5: 20.5, 0.7: 21.5}

for zs in [0.3, 0.5, 0.7, 1.5]:
    # if glob.glob('comb_zs{0}.pkl'.format(zs)) != []:
    #     obs_dict, sim_dict = pickle.load(open(glob.glob('comb_zs{0}.pkl'.format(zs))[0],'rb'))
    # if zs<1:
    #     Imag= Imag_list[zs]
    #     obs_dict, sim_dict = load_HSC_comp(zs =zs, I_mag_break=Imag, no_noise=True)
    #     pickle.dump([obs_dict, sim_dict], open('comb_zs{0}_Imag_{1}_nonoise.pkl'.format(zs,Imag), 'wb'))    
    # elif zs>1: 
    #     obs_dict, sim_dict = load_HST_comp(no_noise=False)
    #     pickle.dump([obs_dict, sim_dict], open('comb_zs{0}.pkl'.format(zs), 'wb'))    
    # print(np.std(sim_dict['Illustris']['BH_Mass_nois_sl']- (m_ml*sim_dict['Illustris']['Stellar_Mass_nois_sl']+b_ml)))
    # plt.scatter(sim_dict['Illustris']['Stellar_Mass'], sim_dict['Illustris']['BH_Mass_nois_sl'])
    if zs<1:
        Imag= Imag_list[zs]
        obs_dict, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}_Imag_{1}.pkl'.format(zs,Imag))[0],'rb'))
    else:
        obs_dict, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}.pkl'.format(zs))[0],'rb'))
        
    if zs == 0.5:
        _, sim_dict_ = pickle.load(open(glob.glob('pickles/comb_zs{0}_Imag_{1}.pkl'.format(zs,21.0))[0],'rb'))
        sim_dict['MBII'] = sim_dict_['MBII']
    
    for sname in ['SAM', 'MBII', 'Illustris', 'TNG100', 'TNG300', 'Horizon']:
        sim = sim_dict[sname]
        imf = imf_dict[sname]
        obs = obs_dict[imf]
        
        if sname == 'MBII':
            if zs==0.5: 
                obs = obs_dict['zs06']   
                zs_ = 0.6
            if zs == 0.7:
                axs[i,j].axes.yaxis.set_ticklabels([])
                axs[i,j].axes.xaxis.set_ticklabels([])
                axs[i,j].xaxis.set_visible(False) 
                axs[i,j].yaxis.set_visible(False) 
                i = i+1
                continue
        else:
            zs_ = zs
        
        sim['BH_Mass'] = sim['BH_Mass'][sim['Stellar_Mass']>0]
        sim['Stellar_Mass'] = sim['Stellar_Mass'][sim['Stellar_Mass']>0]
        
        panel=axs[i,j].hist2d(sim['Stellar_Mass'], sim['BH_Mass'],
                          norm=mpl.colors.LogNorm(), density = True, cmap='summer',bins=50,zorder=-1,
                          alpha=0.5, cmin = 0.001)
        if i == 0:
            sns.kdeplot(sim['Stellar_Mass'], sim['BH_Mass'], linewidths = 2*eps, color = 'seagreen', 
                        fill=True, alpha=0.5, zorder = -10, ax = axs[i,j])

        s, alpha = 150, 0.7
        if sname == 'Horizon':
            sname = 'Horizon-AGN'
        axs[i,j].scatter(sim['Stellar_Mass_nois_sl'][0], sim['BH_Mass_nois_sl'][0],c=colors_sim[i],
                    s=s*eps, marker=".",zorder=0, 
                    edgecolors='black', alpha = alpha, label='{1} z={0}'.format(zs_, sname))
        axs[i,j].scatter(obs['Mstar'][0], #[obs['HSC_ps_mag']<22],
                    obs['MBHs'][0], #[obs['HSC_ps_mag']<22],
                    c='orange',
                    s=s*eps, marker=".",zorder=0.5, edgecolors='black', 
                    alpha = alpha, label=obsname[zs>1])
        #Random place the zorder for zs 0.5 and zs 0.7
        scatters_x = np.concatenate((sim['Stellar_Mass_nois_sl'][1:500], obs['Mstar'][1:])) 
        scatters_y = np.concatenate((sim['BH_Mass_nois_sl'][1:500], obs['MBHs'][1:])) 
        colors = [colors_sim[i]]*len(sim['Stellar_Mass_nois_sl'][1:500]) + ['orange']*len(obs['MBHs'][1:])
        colors = np.asanyarray(colors)
        orders = np.arange(len(scatters_x))
        if zs ==0.5 or zs == 0.7:
            np.random.shuffle(orders)
        axs[i,j].scatter(scatters_x[orders], 
                    scatters_y[orders], 
                    c=colors[orders],
                    s=s*eps, marker=".",zorder=0.5, edgecolors='black', 
                    alpha = alpha)
        xl = np.linspace(5, 13, 100)
        detlaM  = 0
        if imf == 'Sal': 
            detlaM = 0.23
        axs[i,j].plot(xl+detlaM, m_ml*xl+b_ml, color="k", linewidth=1.5*eps,zorder=-0.5)
        # plt.tick_params(which='both', width=2, top=True, right=True,direction='in')
        axs[i,j].grid(linestyle='--')
        axs[i,j].tick_params(labelsize=15*eps)
        axs[i,j].tick_params(which='major', length=7*eps,direction='in')
        axs[i,j].tick_params(which='minor', length=5*eps,direction='in')
        axs[i,j].legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':14*eps},ncol=1, 
                        handletextpad=0.5, handlelength=0, frameon=False )
        if zs<1:
            axs[i,j].set_xlim(9,12.5)
            axs[i,j].set_ylim(6.0,10.3)
        elif zs>1:
            axs[i,j].set_xlim(9.7,11.9)  #
            axs[i,j].set_ylim(7.2, 9.4)  #
        axs[i,j].xaxis.set_minor_locator(AutoMinorLocator())
        axs[i,j].yaxis.set_minor_locator(AutoMinorLocator())
        if i == 5:
            axs[i,j].set_xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=18*eps)
        if j == 0:
            axs[i,j].set_ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=18*eps)
        if i != 5:
            axs[i,j].axes.xaxis.set_ticklabels([])
        if j ==1 or j==2:
            axs[i,j].axes.yaxis.set_ticklabels([])
        if j == 3:
            box = axs[i,3].get_position()
            box.x0 = box.x0 + 0.03
            box.x1 = box.x1 + 0.03
            axs[i,3].set_position(box)
        i = i+1
    j = j+1
    i = 0
plt.savefig('MM_sum.png')
plt.show()    

# # %%Plot Delta M relation
# import copy
# sname_l = ['SAM', 'MBII', 'Illustris', 'TNG100', 'Horizon']
# # cal_M_range = np.arange(9., 12.51, 0.4)
# for zs in [0.7]:
#     # if glob.glob('comb_zs{0}.pkl'.format(zs)) != []:
#     #     obs_dict, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}.pkl'.format(zs))[0],'rb'))
#     # else:
#     #     if zs<1:
#     #         obs_dict, sim_dict = load_HSC_comp(zs =zs)
#     #     elif zs>1: 
#     #         obs_dict, sim_dict = load_HST_comp()
#     #     pickle.dump([obs_dict, sim_dict], open('pickles/comb_zs{0}.pkl'.format(zs), 'wb'))   
    
#     if zs<1:
#         Imag= Imag_list[zs]
#         obs_dict, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}_Imag_{1}.pkl'.format(zs,Imag))[0],'rb'))
#     else:
#         obs_dict, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}.pkl'.format(zs))[0],'rb'))
        
#     if zs == 0.5:
#         _, sim_dict_ = pickle.load(open(glob.glob('pickles/comb_zs{0}_Imag_{1}.pkl'.format(zs,21.0))[0],'rb'))
#         sim_dict['MBII'] = sim_dict_['MBII']
    
#     for idx in range(len(sname_l)):
#     # for idx in [3]:
#         sname = sname_l[idx]
#         obsname = ['HSC', 'HST']
#         colors_sim = ['deepskyblue', 'steelblue', 'c', 'deeppink', 'm']
#         zs_ = zs
#         if sname == 'MBII':
#             if zs==0.5 or zs == 0.7: 
#                 obs = obs_dict['zs06']   
#                 zs_ = 0.6
#         imf = imf_dict[sname]
#         detlaM  = 0
#         if imf == 'Sal': 
#             detlaM = 0.23
#         f,ax=plt.subplots(1,2,figsize=(12,10),gridspec_kw={'width_ratios': [7, 1]}, sharey = True)
#         obj=ax[0]
#         sm_int, bh_int = sim_dict[sname]['Stellar_Mass']- detlaM, sim_dict[sname]['BH_Mass']
#         sm_sim, bh_sim = sim_dict[sname]['Stellar_Mass_nois_sl'][:500]- detlaM, sim_dict[sname]['BH_Mass_nois_sl'][:500]
#         sm_obs, bh_obs = obs_dict[imf]['Mstar'] - detlaM,obs_dict[imf]['MBHs']
        
#         off_int = sm_int, bh_int - (m_ml*sm_int+b_ml)
#         off_sim = sm_sim, bh_sim - (m_ml*sm_sim+b_ml),
#         off_obs = sm_obs, bh_obs - (m_ml*sm_obs+b_ml),
#         panel2=obj.hist2d(off_int[0], off_int[1],
#                           norm=mpl.colors.LogNorm(), density = True, cmap='summer',bins=50,zorder=-1,
#                               alpha=0.5, cmin = 0.001)# , cmax = 1.1)
        
#         s, alpha = 420, 0.7
#         if zs > 1:
#             s,alpha= 800, 0.9
#         # cal_M_range = np.arange(9., 12.51, 0.4)
#         # if zs == 0.3:
#         #     cal_M_range = cal_M_range[1:]
#         #     cal_M_range[0] = 9.5
#         if zs <1 :
#             obs_scatter, sim_scatter, int_scatter = [], [], []
            
#             sm_ = copy.deepcopy(sm_obs)
#             sm_.sort()
#             m_min, m_max = np.mean(sm_[:3]), np.mean(sm_[-3:])
#             cal_M_range_0 = np.arange(m_min, m_max, 0.3)
#             for i in range(len(cal_M_range_0)-1):
#                 s_bool = (sm_obs>cal_M_range_0[i])*(sm_obs<cal_M_range_0[i+1])
#                 cal_HSC_Mstar = sm_obs[s_bool]
#                 cal_HSC_MBHs = bh_obs[s_bool]
#                 obs_res = cal_HSC_MBHs-(m_ml*cal_HSC_Mstar+b_ml)
#                 obs_scatter.append( [np.mean(obs_res), np.std(obs_res)] )
#             sm_ = copy.deepcopy(sm_sim)
#             sm_.sort()
#             m_min, m_max = np.mean(sm_[:3]), np.mean(sm_[-10:])
#             cal_M_range_1 = np.arange(m_min, m_max, 0.3)
#             for i in range(len(cal_M_range_1)-1):
#                 s_bool = (sm_sim>cal_M_range_1[i])*(sm_sim<cal_M_range_1[i+1])
#                 cal_Mstar = sm_sim[s_bool]
#                 cal_MBHs = bh_sim[s_bool]
#                 obs_res = cal_MBHs-(m_ml*cal_Mstar+b_ml)
#                 sim_scatter.append( [np.mean(obs_res), np.std(obs_res)] )
#             sm_ = copy.deepcopy(sm_int)
#             sm_.sort()
#             m_min, m_max = np.mean(sm_[:3]), np.mean(sm_[-3:])
#             if sname == 'SAM' and zs == 0.3:
#                 m_min = 9.3
#             cal_M_range_2 = np.arange(m_min, m_max, 0.3)
#             for i in range(len(cal_M_range_2)-1):
#                 s_bool = (sm_int>cal_M_range_2[i])*(sm_int<cal_M_range_2[i+1])
#                 cal_Mstar = sm_int[s_bool]
#                 cal_MBHs = bh_int[s_bool]
#                 obs_res = cal_MBHs-(m_ml*cal_Mstar+b_ml)
#                 int_scatter.append( [np.mean(obs_res), np.std(obs_res)] )    
#             obs_scatter = np.array(obs_scatter)
#             sim_scatter = np.array(sim_scatter)
#             int_scatter = np.array(int_scatter)
#             ax[0].errorbar(cal_M_range_2[:-1]+ (cal_M_range_2[1]-cal_M_range_2[0])/2+0.15, int_scatter[:,0], int_scatter[:,1], color = 'black', 
#                   zorder = 10, linewidth = 3, linestyle= '-',fmt='o', alpha = 0.9)
#             ax[0].scatter(cal_M_range_2[:-1]+ (cal_M_range_2[1]-cal_M_range_2[0])/2+0.15, int_scatter[:,0], s = 300, color = 'darkseagreen', marker = 's',
#                           edgecolor = 'black', linewidth=3, zorder = 11)    
#             ax[0].errorbar(cal_M_range_0[:-1]+ (cal_M_range_0[1]-cal_M_range_0[0])/2, obs_scatter[:,0], obs_scatter[:,1], color = 'black', 
#                   zorder = 10, linewidth = 3, linestyle= '-',fmt='o', alpha = 0.9)
#             ax[0].scatter(cal_M_range_0[:-1]+ (cal_M_range_0[1]-cal_M_range_0[0])/2, obs_scatter[:,0], s = 300, color = 'orange', marker = 's',
#                           edgecolor = 'black', linewidth=3, zorder = 11)
#             ax[0].errorbar(cal_M_range_1[:-1]+ (cal_M_range_1[1]-cal_M_range_1[0])/2+0.05, sim_scatter[:,0], sim_scatter[:,1], color = 'black', 
#                   zorder = 10, linewidth = 3, linestyle= '-',fmt='o', alpha = 0.9)
#             ax[0].scatter(cal_M_range_1[:-1]+ (cal_M_range_1[1]-cal_M_range_1[0])/2+0.05, sim_scatter[:,0], s = 300, color = colors_sim[idx], marker = 's',
#                           edgecolor = 'black', linewidth=3, zorder = 11)    
#         ax[0].plot(np.linspace(7, 13, 100), np.linspace(7, 13, 100) *0, 'k', zorder = 100, linewidth = 3 )
#         ssname = sname
#         if sname == 'Horizon':
#             ssname = 'Horizon-AGN'
#         ax[0].scatter(off_sim[0], off_sim[1],
#                     c=colors_sim[idx],
#                     s=420, marker=".",zorder=0, edgecolors='k', alpha = 0.7, label='{1} sample z={0}'.format(zs_, ssname))
#         ax[0].scatter(off_obs[0], off_obs[1],
#                     c='orange',
#                     s=s, marker=".",zorder=1, edgecolors='k', alpha = alpha, label=obsname[zs>1])
        
#         ax[0].set_xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
#         ax[0].set_ylabel(r"$\Delta$logM$_{\rm BH}$ (vs M$_*$)",fontsize=35)
        
#         if zs<1:
#             ax[0].set_xlim(9,12.5)
#             ax[0].set_ylim(-1.7, 2.5)
#         elif zs>1:
#             ax[0].set_xlim(9.7,11.5)  #
#             ax[0].set_ylim(-1.0, 1.7)  #
#         ax[0].grid(linestyle='--')
#         ax[0].tick_params(labelsize=25)
#         ax[0].tick_params(which='both', width=2, top=True, right=True,direction='in')
#         ax[0].tick_params(which='major', length=10)
#         ax[0].tick_params(which='minor', length=6)#, color='r’)
#         ax[0].legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':32},ncol=2,#handletextpad=0,
#                       handletextpad=0.5, handlelength=0,)
#         ax[0].xaxis.set_minor_locator(AutoMinorLocator())
#         ax[0].yaxis.set_minor_locator(AutoMinorLocator())
#         his_xy0_ =  ax[1].hist(off_int[1], orientation='horizontal'
#                     , histtype=u'step',density=True, color = 'darkseagreen', linewidth = 4)
#         his_xy1_ = ax[1].hist(off_sim[1], orientation='horizontal'
#                     , histtype=u'step',density=True, color = colors_sim[idx], linewidth = 4)
#         his_xy2_ = ax[1].hist(off_obs[1], orientation='horizontal'
#                     , histtype=u'step',density=True, color = 'orange', linewidth = 4)
#         his_max = np.max([his_xy0_[0].max(), his_xy1_[0].max(), his_xy2_[0].max()])
#         sim_mean = np.mean(off_sim[1])
#         obs_mean = np.mean(off_obs[1])
#         int_mean = np.mean(off_int[1])
#         ax[1].plot([0, 10], [int_mean,int_mean], ls = '--',  linewidth = 3,color = 'darkseagreen', zorder = 0)
#         ax[1].plot([0, 10], [sim_mean,sim_mean], linewidth = 3,color = colors_sim[idx], zorder = 1)
#         ax[1].plot([0, 10], [obs_mean,obs_mean], linewidth = 3,color = 'orange', zorder = 0.5)
#         ax[1].set_xlim(0, his_max*1.2)
#         ax[1].set_xticks([])
        
#         f.tight_layout()
#         plt.subplots_adjust(wspace=0.01)
#        # plt.savefig('DeltaMM_{1}_zs_{0}.png'.format(str(zs_).replace('.',''), sname))
#         plt.show()
            
