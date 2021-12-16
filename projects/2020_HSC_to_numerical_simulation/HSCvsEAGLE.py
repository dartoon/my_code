#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 15:21:57 2021

@author: Dartoon
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import scipy.stats as st
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

#%%
from prep_comparison import HSC_set, EAGLE_set, comp_plot

# filenames = glob.glob('EAGLE/*') 
# filenames.sort()
zs = 0.7

filename = 'EAGLE/simdata_{0}.txt'.format(zs)

# filename = filenames[3]
# zs = 0.5

# filename = filenames[4]
# zs = 0.3    #Actually from z = 0, not work, too few sample

if zs < 0.5:
    # HSC_Mstar = HSC_Mstar_overall[HSC_z<0.5]
    # HSC_MBHs = HSC_MBHs_overall[HSC_z<0.5]
    I_mag_break = 21  #z~0.3 #!!!
if zs >= 0.5:    
    # HSC_Mstar = HSC_Mstar_overall[HSC_z>0.5]
    # HSC_MBHs = HSC_MBHs_overall[HSC_z>0.5]
    I_mag_break = 22.0    #z~0.7

detlaM  = 0
imf = 'Sal'
if imf == 'Sal':
    detlaM = 0.23

ifplot = True

    #%%

for ii in range(1):
    HSC = HSC_set(zs, core = True,imf = imf)
    EAGLE = EAGLE_set(filename, HSC_Lbol_overall=HSC['HSC_Lbol_overall'], HSC_MBHs_overall=HSC['HSC_MBHs_overall'],
                  zs = zs, I_mag_break = I_mag_break, imf =  imf)
    
    #%%
    comp_plot(EAGLE['BH_Mass'], EAGLE['Eddington_ratio'], 'BH_Mass', 'Eddington_ratio')
    # plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
    # plt.ylabel('AGN Eddington ratio', fontsize=25)
    plt.tick_params(labelsize=25)
    plt.xlim(5,10)
    if ifplot == True:
        plt.show()
    else:
        plt.close()
        
        
    comp_plot(EAGLE['logLbol_nois'], EAGLE['BH_Mass_nois'], 'logLbol_nois', 'BH_Mass_nois', alpha = 0.2, label = 'EAGLE')
    plt.scatter(EAGLE['logLbol_nois_sl'], EAGLE['BH_Mass_nois_sl'], color = 'green',alpha=0.2, zorder = 1, label = 'selected EAGLE')
    plt.scatter(HSC['HSC_Lbol_overall'], HSC['HSC_MBHs_overall'],c='orange',alpha=0.2,zorder = 0.5, label = 'HSC QSO distribution')
    plt.xlabel('AGN logLbol', fontsize=25)
    plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
    plt.tick_params(labelsize=25)
    plt.xlim(30, 46.5)
    plt.ylim(5.8,10)
    plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':20},ncol=2,handletextpad=0)
    if ifplot == True:
        plt.show()
    else:
        plt.close()
        
    #%%
    # plt.figure(figsize=(11.5,12))      
    # plt.scatter(EAGLE['Stellar_Mass_nois'], EAGLE['BH_Mass_nois'],c='gray',
    #             s=220, marker=".",zorder=-10, edgecolors='k', alpha = 0.2)
    # plt.scatter(EAGLE['Stellar_Mass_nois_sl'], EAGLE['BH_Mass_nois_sl'],c='m',
    #             s=220, marker=".",zorder=0, edgecolors='k', alpha = 0.7, label='EAGLE-AGN sample z={0}'.format(zs))
    # # plt.scatter(HSC['HSC_Mstar'],HSC['HSC_MBHs'],c='orange',
    # #             s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')
    # plt.scatter(HSC['HSC_Mstar'][HSC['HSC_ps_mag']<I_mag_break],HSC['HSC_MBHs'][HSC['HSC_ps_mag']<I_mag_break],c='orange',
    #             s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')
    
    # m_ml, b_ml = (0.981139684856507, -2.545890295477823)
    # xl = np.linspace(5, 13, 100)
    # plt.plot(xl+ detlaM, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)
    # plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
    # plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
    # plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
    # plt.xlim(9,12.5)
    # plt.ylim(6.0,10.3)
    # plt.grid(linestyle='--')
    # plt.tick_params(labelsize=25)
    # plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
    # if ifplot == True:
    #     plt.show()
    # else:
    #     plt.close()
        
    f,ax=plt.subplots(1,1,figsize=(14,12))   
    import matplotlib as mpl
    obj=ax
    

    s, alpha = 420, 0.7
    if zs == 0.3:
        s = 620
        alpha = 0.9   
    
    # panel2=obj.hist2d(EAGLE['Stellar_Mass_reali'], EAGLE['BH_Mass_reali'], #!!!
    panel2=obj.hist2d(EAGLE['Stellar_Mass'], EAGLE['BH_Mass'],
                      norm=mpl.colors.LogNorm(), density = True, cmap='summer',bins=50,zorder=-1,
                      alpha=0.5, cmin = 0.001)# , cmax = 1.1)
    
    plt.scatter(EAGLE['Stellar_Mass_nois_sl'][:500], EAGLE['BH_Mass_nois_sl'][:500],c='m',
                s=420, marker=".",zorder=0, edgecolors='k', alpha = 0.7, label='EAGLE-AGN sample z={0}'.format(zs))
    # plt.scatter(HSC['HSC_Mstar'],HSC['HSC_MBHs'],c='orange',
    #             s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')
    plt.scatter(HSC['HSC_Mstar'][HSC['HSC_ps_mag']<I_mag_break][:500],HSC['HSC_MBHs'][HSC['HSC_ps_mag']<I_mag_break][:500],c='orange',
                s=s, marker=".",zorder=-1, edgecolors='k', alpha = alpha, label='HSC sample')
    
    xl = np.linspace(5, 13, 100)
    m_ml, b_ml = (0.981139684856507, -2.545890295477823)
    plt.plot(xl+detlaM, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)
    # plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
    plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
    plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
    plt.xlim(9,12.5)
    plt.ylim(6.0,10.3)
    plt.grid(linestyle='--')
    plt.tick_params(labelsize=25)
    plt.tick_params(which='both', width=2, top=True, right=True,direction='in')
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=6)#, color='r’)
    plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':32},ncol=1,handletextpad=0)
    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    cbar=f.colorbar(panel2[3],ax=obj)
    cbar.ax.tick_params(labelsize=30) 
    plt.savefig('MM_EAGLE_zs_{0}.png'.format(zs))
    if ifplot == True:
        plt.show()
    else:
        plt.close()   

    #%%Plot offset        
    f,ax=plt.subplots(1,2,figsize=(12,10),gridspec_kw={'width_ratios': [7, 1]}, sharey = True)
    # f.suptitle(r"Offset VS M*, z={0}".format(zs), fontsize = 20)
    obj=ax[0]
    
    sm_int, bh_int = EAGLE['Stellar_Mass']- detlaM, EAGLE['BH_Mass']
    sm_sim, bh_sim = EAGLE['Stellar_Mass_nois_sl']- detlaM, EAGLE['BH_Mass_nois_sl']
    sm_obs, bh_obs = HSC['HSC_Mstar'][HSC['HSC_ps_mag']<I_mag_break]- detlaM, HSC['HSC_MBHs'][HSC['HSC_ps_mag']<I_mag_break]
    
    off_int = sm_int, bh_int - (m_ml*sm_int+b_ml)
    off_sim = sm_sim, bh_sim - (m_ml*sm_sim+b_ml),
    off_obs = sm_obs, bh_obs - (m_ml*sm_obs+b_ml),
    panel2=obj.hist2d(off_int[0], off_int[1],
                      norm=mpl.colors.LogNorm(), density = True, cmap='summer',bins=50,zorder=-1,
                          alpha=0.5, cmin = 0.001)# , cmax = 1.1)

    s, alpha = 420, 0.7
    cal_M_range = np.arange(9.5, 12.1, 0.3)
    if zs == 0.3:
        cal_M_range = cal_M_range[1:]
        cal_M_range[0] = 9.5
        s = 620
        alpha = 0.9
    obs_scatter, sim_scatter = [], []
    for i in range(len(cal_M_range)-1):
        s_bool = (sm_obs>cal_M_range[i])*(sm_obs<cal_M_range[i+1])
        cal_HSC_Mstar = sm_obs[s_bool]
        cal_HSC_MBHs = bh_obs[s_bool]
        obs_res = cal_HSC_MBHs-(m_ml*cal_HSC_Mstar+b_ml)
        obs_scatter.append( [np.mean(obs_res), np.std(obs_res)] )
        s_bool = (sm_sim>cal_M_range[i])*(sm_sim<cal_M_range[i+1])
        cal_HSC_Mstar = sm_sim[s_bool]
        cal_HSC_MBHs = bh_sim[s_bool]
        obs_res = cal_HSC_MBHs-(m_ml*cal_HSC_Mstar+b_ml)
        sim_scatter.append( [np.mean(obs_res), np.std(obs_res)] )
    obs_scatter = np.array(obs_scatter)
    sim_scatter = np.array(sim_scatter)
    ax[0].errorbar(cal_M_range[:-1]+ (cal_M_range[1]-cal_M_range[0])/2, obs_scatter[:,0], obs_scatter[:,1], color = 'black', 
          zorder = 10, linewidth = 3, linestyle= '-',fmt='o', alpha = 0.9)
    ax[0].scatter(cal_M_range[:-1]+ (cal_M_range[1]-cal_M_range[0])/2, obs_scatter[:,0], s = 300, color = 'orange', marker = 's',
                  edgecolor = 'black', linewidth=3, zorder = 11)
    ax[0].errorbar(cal_M_range[:-1]+ (cal_M_range[1]-cal_M_range[0])/2+0.05, sim_scatter[:,0], sim_scatter[:,1], color = 'black', 
          zorder = 10, linewidth = 3, linestyle= '-',fmt='o', alpha = 0.9)
    ax[0].scatter(cal_M_range[:-1]+ (cal_M_range[1]-cal_M_range[0])/2+0.05, sim_scatter[:,0], s = 300, color = 'm', marker = 's',
                  edgecolor = 'black', linewidth=3, zorder = 11)    
    ax[0].plot(np.linspace(7, 13, 100), np.linspace(7, 13, 100) *0, 'k', zorder = 1, linewidth = 3 )
    
    ax[0].scatter(off_sim[0], off_sim[1],
                c='m',
                s=420, marker=".",zorder=0, edgecolors='k', alpha = 0.7, label='EAGLE-AGN sample z={0}'.format(zs))
    ax[0].scatter(off_obs[0], off_obs[1],
                c='orange',
                s=420, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')
    # xl = np.linspace(5, 13, 100)
    # plt.plot(xl, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)
    # plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
    ax[0].set_xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
    ax[0].set_ylabel(r"$\Delta$logM$_{\rm BH}$ (vs M$_*$)",fontsize=35)
    ax[0].set_xlim(9,12.5)
    ax[0].set_ylim(-1.5, 2.4)
    ax[0].grid(linestyle='--')
    ax[0].tick_params(labelsize=25)
    ax[0].tick_params(which='both', width=2, top=True, right=True,direction='in')
    ax[0].tick_params(which='major', length=10)
    ax[0].tick_params(which='minor', length=6)#, color='r’)
    ax[0].legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':32},ncol=1,handletextpad=0)
    ax[0].xaxis.set_minor_locator(AutoMinorLocator())
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    
    his_xy0_ =  ax[1].hist(off_int[1], orientation='horizontal'
               , histtype=u'step',density=True, color = 'green', linewidth = 4)
    his_xy1_ = ax[1].hist(off_sim[1], orientation='horizontal'
               , histtype=u'step',density=True, color = 'm', linewidth = 4)
    his_xy2_ = ax[1].hist(off_obs[1], orientation='horizontal'
               , histtype=u'step',density=True, color = 'orange', linewidth = 4)
    his_max = np.max([his_xy0_[0].max(), his_xy1_[0].max(), his_xy2_[0].max()])
    # ax[1].set_yticks([])
    sim_mean = np.mean(off_sim[1])
    obs_mean = np.mean(off_obs[1])
    ax[1].plot([0, 10], [sim_mean,sim_mean], linewidth = 3,color = 'm')
    ax[1].plot([0, 10], [obs_mean,obs_mean], linewidth = 3,color = 'orange')
    ax[1].set_xlim(0, his_max*1.2)
    ax[1].set_xticks([])
    
    f.tight_layout()
    plt.subplots_adjust(wspace=0.01)
    # from matplotlib.ticker import AutoMinorLocator
    # cbar=f.colorbar(panel2[3],ax=obj)
    # cbar.ax.tick_params(labelsize=30) 
    plt.savefig('DeltaMM_EAGLE_zs_{0}.png'.format(zs))
    if ifplot == True:
        plt.show()
    else:
        plt.close()
    
    cals = off_int[1]#[(off_int[0]<off_obs[0].max())*(off_int[0]>off_obs[0].min())]
    print('{0:.2f}, {1:.2f}'.format(np.mean(cals), np.std(cals)))
            
    #%%
    # EAGLE_scatter = (EAGLE['BH_Mass_nois_sl'] - ( m_ml* (EAGLE['Stellar_Mass_nois_sl'] - detlaM)+b_ml ) )
    # EAGLE_scatter_nosl = (EAGLE['BH_Mass_nois'] - ( m_ml* (EAGLE['Stellar_Mass_nois'] -detlaM) +b_ml ) )
    # HSC_scatter = (HSC['HSC_MBHs'] - ( m_ml*(HSC['HSC_Mstar'] - detlaM)+b_ml ) )
    # #Plot the 1-D scatter for MM.
    # fig, ax = plt.subplots(figsize=(8,7))
    # plt.hist(EAGLE_scatter_nosl, histtype=u'step',density=True,
    #           label=('EAGLE-AGN sample scatter nosl'), linewidth = 2, color='gray')
    # plt.hist(EAGLE_scatter,histtype=u'step',density=True,
    #           label=('EAGLE-AGN sample scatter'), linewidth = 2, color='m')
    # plt.hist(HSC_scatter, histtype=u'step',density=True,
    #           label=('HSC sample scatter'), linewidth = 2, color='orange')
    # # plt.hist(EAGLE_scatter_noselect,histtype=u'step',density=True,
    # #           label=('EAGLE-AGN sample scatter, no selection'), linewidth = 2, color='gray')
    # plt.title(r"The offset comparison for the M$_{\rm BH}$-M$_{*}$ relation", fontsize = 25)
    # plt.tick_params(labelsize=20)
    # plt.legend(prop={'size':10})
    # plt.yticks([])
    # # ax.xaxis.set_minor_locator(AutoMinorLocator())
    # plt.tick_params(which='both', width=2, top=True,direction='in')
    # plt.tick_params(which='major', length=10)
    # plt.tick_params(which='minor', length=6)#, color='r’)
    # plt.xlabel(r'$\Delta$log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
    # #plt.savefig('comp_scatter_MM_EAGLEonly.pdf')
    # if ifplot == True:
    #     plt.show()
    # else:
    #     plt.close()

    # sim_offset_nosl = EAGLE_scatter_nosl 
    # sim_offset = EAGLE_scatter
    # obs_offset = HSC_scatter
    # rfilename = 'offset_result/' + 'EAGLE_zs{0}.txt'.format(zs)
    # if_file = glob.glob(rfilename)
    # if ii == 0:
    #     write_file =  open(rfilename,'w') 
    # for i in range(max(len(sim_offset), len(obs_offset))):
    #     try:
    #         write_file.write('{0} {1} {2} {3} {4} {5} {6}'.format(sim_offset_nosl[i], sim_offset[i], obs_offset[i], 
    #                                                               EAGLE['Stellar_Mass_nois_sl'][i], EAGLE['BH_Mass_nois_sl'][i], HSC['HSC_Mstar'][i], HSC['HSC_MBHs'][i] ))
    #     except:
    #         try:
    #             write_file.write('{0} {1} -99 {2} {3} -99 -99'.format(sim_offset_nosl[i], sim_offset[i], EAGLE['Stellar_Mass_nois_sl'][i], EAGLE['BH_Mass_nois_sl'][i]))
    #         except:            
    #             write_file.write('{0} -99 {1} -99 -99 {2} {3}'.format(sim_offset_nosl[i], obs_offset[i],HSC['HSC_Mstar'][i], HSC['HSC_MBHs'][i]))
    #     write_file.write("\n")
    # # from scipy import stats
    # sim_scatter_std = np.std(EAGLE_scatter)
    # obs_scatter_std = np.std(HSC_scatter)
    # if ii%50 == 0:
    #     print(ii)    
# write_file.close()