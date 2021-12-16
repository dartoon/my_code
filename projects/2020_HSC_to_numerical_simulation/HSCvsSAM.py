#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 10:27:56 2021

@author: Dartoon
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

m_ml, b_ml = (0.981139684856507, -2.545890295477823)

from prep_comparison import SAM_set, HSC_set, comp_plot, quasar_filter

ifplot = True

zs = 0.3
imf = 'Sal'
HSC = HSC_set(zs, core=True, imf = imf)
# HSC['HSC_Mstar_overall'] = HSC['HSC_Mstar_overall']+ 0.25
# HSC['HSC_Mstar'] = HSC['HSC_Mstar']+ 0.25
if zs == 0.7:
    filename = 'SAM/catalogue2.dat'
else:
    filename = 'SAM/catalogue.dat'

if zs <= 0.5:
    # HSC_Mstar = HSC_Mstar_overall[HSC_z<0.5]
    # HSC_MBHs = HSC_MBHs_overall[HSC_z<0.5]
    I_mag_break = 20.5  #z~0.3
if zs > 0.5:    
    # HSC_Mstar = HSC_Mstar_overall[HSC_z>0.5]
    # HSC_MBHs = HSC_MBHs_overall[HSC_z>0.5]
    I_mag_break = 22.0    #z~0.7
    
    
for ii in range(1):    
    SAM = SAM_set(filename, zs=zs, HSC_Lbol_overall=HSC['HSC_Lbol_overall'], HSC_MBHs_overall=HSC['HSC_MBHs_overall'],
                  I_mag_break = [20.5,22.0 ])
    
    detlaM  = 0
    if imf == 'Sal':
        detlaM = 0.23
    
    SAM_scatter_overall = (SAM['BH_Mass_nois'] - ( m_ml* ( SAM['Stellar_Mass_nois'] - detlaM) +b_ml ) )
    SAM_scatter = (SAM['BH_Mass_nois_sl'] - ( m_ml* ( SAM['Stellar_Mass_nois_sl'] - detlaM) +b_ml ) )
    
    # text  = np.loadtxt(filename)
    # text = text[text[:,5] != 0]
    # redshift = text[:, 0]
    # logMBH =  text[:, 1]
    # logM_mass =  text[:, 2]
    # g_galaxy =  text[:, 3]
    # r_galaxy =  text[:, 4]
    # AGN_bol_10_45 =  text[:, 5]
    # stat_weight =  text[:, 6]
    
    # plt.scatter(redshift, stat_weight)
    # plt.xlabel("$z$", fontsize = 15)
    # plt.ylabel("$stat_weight$", fontsize = 15)
    # plt.show()
    
    # plt.scatter(redshift, AGN_bol_10_45)
    # plt.xlabel("$z$", fontsize = 15)
    # plt.ylabel("Lbol (10^45)", fontsize = 15)
    # plt.show()
    
    # plt.scatter(redshift, logMBH)
    # plt.xlabel("$z$", fontsize = 15)
    # plt.ylabel("$logMBH$", fontsize = 15)
    # plt.show()
    
    
    #%%
    # plt.figure(figsize=(12.5,12))      
    # plt.scatter(SAM['Stellar_Mass_nois'], SAM['BH_Mass_nois'],c='gray',
    #             s=220, marker=".",zorder=1, edgecolors='k', alpha = 0.2, label='SAM sample no sl')
    
    # plt.scatter(SAM['Stellar_Mass_nois_sl'][:len(HSC['HSC_Mstar'])], SAM['BH_Mass_nois_sl'][:len(HSC['HSC_Mstar'])],c='green',
    #             s=220, marker=".",zorder=1, edgecolors='k', alpha = 0.7, cmap='autumn', label='SAM sample z={0}'.format(zs))
    
    # plt.scatter(HSC['HSC_Mstar'],HSC['HSC_MBHs'],c='orange',
    #             s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')
    
    # xl = np.linspace(5, 13, 100)
    
    
    # plt.plot(xl + detlaM, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)
    # plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
    # plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
    # plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
    # plt.xlim(9,12.5)
    # plt.ylim(6.0,10.3)
    # plt.grid(linestyle='--')
    # plt.tick_params(labelsize=25)
    # plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=1,handletextpad=0)
    # plt.savefig('MM_SAM_zs_{0}.pdf'.format(zs))
    # if ifplot == True:
    #     plt.show()
    # else:
    #     plt.close()
    
    f,ax=plt.subplots(1,1,figsize=(14,12))   
    import matplotlib as mpl
    obj=ax
    # panel2=obj.hist2d(SAM['Stellar_Mass_reali'], SAM['BH_Mass_reali'], #!!!
    panel2=obj.hist2d(SAM['Stellar_Mass_reali'], SAM['BH_Mass_reali'],
                      norm=mpl.colors.LogNorm(), density = True, cmap='summer',bins=50,zorder=10,
                      alpha=0.5, cmin = 0.001 , cmax = 1.1)
    
    # import seaborn as sns
    # sns.kdeplot(SAM['Stellar_Mass_reali'], SAM['BH_Mass_reali'], linewidths = 2, color = 'green', fill=True, alpha=0.5, zorder = -10)

    s, alpha = 420, 0.7
    cal_M_range = np.arange(9.5, 12.1, 0.3)
    if zs == 0.3:
        cal_M_range = cal_M_range[1:]
        cal_M_range[0] = 9.5
        s = 620
        alpha = 0.9
    
    plt.scatter(SAM['Stellar_Mass_nois_sl'][:500], SAM['BH_Mass_nois_sl'][:500],c='pink',
                s=420, marker=".",zorder=1.2, edgecolors='k', alpha = 0.7, label='SAM sample z={0}'.format(zs))
    # plt.scatter(HSC['HSC_Mstar'],HSC['HSC_MBHs'],c='orange',
    #             s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')
    plt.scatter(HSC['HSC_Mstar'][HSC['HSC_ps_mag']<I_mag_break][:500],HSC['HSC_MBHs'][HSC['HSC_ps_mag']<I_mag_break][:500],c='orange',
                s=s, marker=".",zorder=1.3, edgecolors='k', alpha = alpha, label='HSC sample')
    
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
    cbar=f.colorbar(panel2[3],ax=obj, ticks=[])
    cbar.ax.tick_params(labelsize=30) 
    plt.savefig('MM_SAM_zs_{0}.png'.format(zs))
    plt.show()        
    
    #%%    
    import matplotlib as mpl
    from matplotlib.ticker import AutoMinorLocator
    f,ax=plt.subplots(1,2,figsize=(12,10),gridspec_kw={'width_ratios': [7, 1]}, sharey = True)
    # f.suptitle(r"Offset VS M*, z={0}".format(zs), fontsize = 20)
    obj=ax[0]
    
    sm_int, bh_int = SAM['Stellar_Mass_reali']- detlaM, SAM['BH_Mass_reali']
    sm_sim, bh_sim = SAM['Stellar_Mass_nois_sl'][:500]- detlaM, SAM['BH_Mass_nois_sl'][:500]
    sm_obs, bh_obs = HSC['HSC_Mstar'][HSC['HSC_ps_mag']<I_mag_break][:500]- detlaM,HSC['HSC_MBHs'][HSC['HSC_ps_mag']<I_mag_break][:500]
    
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
    ax[0].scatter(cal_M_range[:-1]+ (cal_M_range[1]-cal_M_range[0])/2+0.05, sim_scatter[:,0], s = 300, color = 'pink', marker = 's',
                  edgecolor = 'black', linewidth=3, zorder = 11)    
    ax[0].plot(np.linspace(7, 13, 100), np.linspace(7, 13, 100) *0, 'k', zorder = 1, linewidth = 3 )
    
    
    ax[0].scatter(off_sim[0], off_sim[1],
                c='pink',
                s=420, marker=".",zorder=0, edgecolors='k', alpha = 0.7, label='SAM sample z={0}'.format(zs))
    ax[0].scatter(off_obs[0], off_obs[1],
                c='orange',
                s=s, marker=".",zorder=1, edgecolors='k', alpha = alpha, label='HSC sample')
    
    # import seaborn as sns
    # sns.kdeplot(off_int[0], off_int[1], linewidths = 1, color = 'green', ax=ax[0])
    
    # xl = np.linspace(5, 13, 100)
    # plt.plot(xl, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)
    # plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
    ax[0].set_xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
    ax[0].set_ylabel(r"$\Delta$logM$_{\rm BH}$ (vs M$_*$)",fontsize=35)
    ax[0].set_xlim(9,12.5)
    ax[0].set_ylim(-1.7, 2.5)
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
               , histtype=u'step',density=True, color = 'pink', linewidth = 4)
    his_xy2_ = ax[1].hist(off_obs[1], orientation='horizontal'
               , histtype=u'step',density=True, color = 'orange', linewidth = 4)
    his_max = np.max([his_xy0_[0].max(), his_xy1_[0].max(), his_xy2_[0].max()])
    # ax[1].set_yticks([])
    sim_mean = np.mean(off_sim[1])
    obs_mean = np.mean(off_obs[1])
    ax[1].plot([0, 10], [sim_mean,sim_mean], linewidth = 3,color = 'pink')
    ax[1].plot([0, 10], [obs_mean,obs_mean], linewidth = 3,color = 'orange')
    ax[1].set_xlim(0, his_max*1.2)
    ax[1].set_xticks([])
    
    f.tight_layout()
    plt.subplots_adjust(wspace=0.01)
    from matplotlib.ticker import AutoMinorLocator
    # cbar=f.colorbar(panel2[3],ax=obj)
    # cbar.ax.tick_params(labelsize=30) 
    plt.savefig('DeltaMM_SAM_zs_{0}.png'.format(zs))
    plt.show()
    
    #%%
    
    cals = off_int[1]#[(off_int[0]<off_obs[0].max())*(off_int[0]>off_obs[0].min())]
    print('{0:.2f}, {1:.2f}'.format(np.mean(cals), np.std(cals)))    
    
    
    # #%%
    HSC_scatter = (HSC['HSC_MBHs'] - ( m_ml*(HSC['HSC_Mstar']-detlaM)+b_ml ) )
    
    # #Plot the 1-D scatter for MM.
    # fig, ax = plt.subplots(figsize=(8,7))
    # plt.hist(SAM_scatter_overall, histtype=u'step',density=True,
    #           label=('SAM sample scatter nosl'), linewidth = 2, color='gray')
    # plt.hist(HSC_scatter, histtype=u'step',density=True,
    #           label=('HSC sample scatter'), linewidth = 2, color='orange')
    # plt.hist(SAM_scatter,histtype=u'step',density=True,
    #           label=('SAM sample scatter'), linewidth = 2, color='green')
    # plt.title(r"The offset comparison for the M$_{\rm BH}$-M$_{*}$ relation", fontsize = 20)
    # plt.tick_params(labelsize=20)
    # plt.legend(prop={'size':10})
    # plt.yticks([])
    
    # # ax.xaxis.set_minor_locator(AutoMinorLocator())
    # plt.tick_params(which='both', width=2, top=True,direction='in')
    # plt.tick_params(which='major', length=10)
    # plt.tick_params(which='minor', length=6)#, color='r’)
    
    # plt.xlabel(r'$\Delta$log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
    # #plt.savefig('comp_scatter_MM_SAMonly.pdf')
    # if ifplot == True:
    #     plt.show()
    # else:
    #     plt.close()
        

    obs_offset = HSC_scatter
    # sim_offset_nosl = SAM_scatter_overall 
    # sim_offset = SAM_scatter
    # rfilename = 'offset_result/' + 'SAM_zs{0}.txt'.format(zs)
    # if_file = glob.glob(rfilename)
    # write_file =  open(rfilename,'w') 
    
    # for i in range(max(len(sim_offset), len(obs_offset))):
    #     try:
    #         write_file.write('{0} {1} {2} {3} {4} {5} {6}'.format(sim_offset_nosl[i], sim_offset[i], obs_offset[i], 
    #                                                               SAM['Stellar_Mass_nois_sl'][i], SAM['BH_Mass_nois_sl'][i], HSC['HSC_Mstar'][i], HSC['HSC_MBHs'][i] ))
    #     except:
    #         try:
    #             write_file.write('{0} {1} -99 {2} {3} -99 -99'.format(sim_offset_nosl[i], sim_offset[i], SAM['Stellar_Mass_nois_sl'][i], SAM['BH_Mass_nois_sl'][i]))
    #         except:            
    #             write_file.write('{0} -99 {1} -99 -99 {2} {3}'.format(sim_offset_nosl[i], obs_offset[i],HSC['HSC_Mstar'][i], HSC['HSC_MBHs'][i]))
    #     write_file.write("\n")
    # write_file.close()
        
    # # from scipy import stats
    # sim_scatter_std = np.std(SAM_scatter)
    # obs_scatter_std = np.std(HSC_scatter)
    # # print("obs scatter:", obs_scatter_std)
    # # print("sim scatter:", sim_scatter_std)
    # # print("KS p-value:", stats.ks_2samp(SAM_scatter, HSC_scatter).pvalue)
    
    # # print("({0:.2f}, {1:.2f})".format( np.mean(SAM_scatter) - np.mean(HSC_scatter), np.std(SAM_scatter) - np.std(HSC_scatter) ))
    
    # # print("for paper Observation", 'zs=', zs)
    print('{0:.2f}, {1:.2f}'.format(np.mean(HSC_scatter), np.std(HSC_scatter)))
    
    # # print("for paper SAM", 'zs=', zs)
    # # print('{0:.2f}, {1:.2f}'.format(np.mean(SAM_scatter), np.std(SAM_scatter)))

    # # rfilename = 'MC_result/' + 'SAM_zs{0}.txt'.format(zs)
    # # if_file = glob.glob(rfilename)
    # # if if_file == []:
    # #     write_file =  open(rfilename,'w') 
    # # else:
    # #     write_file =  open(rfilename,'r+') 
    # #     write_file.read()   
        
    # # write_file.write( "{0:.3f} {1:.3f}".format(np.mean(SAM_scatter), np.std(SAM_scatter)))
    # # write_file.write("\n")
    # # write_file.close()
    # if ii%50 == 0:
    #     print(ii)
        
    
    # #%%
    # comp_plot(HSC['HSC_Mstar'],HSC['HSC_galaxy_abs_iMags'],c = 'orange', alpha=0.2)
    # plt.scatter(SAM['Stellar_Mass_nois_sl'], SAM['sdss_mag_galaxy_sl'], c = 'green',alpha=0.2)
    # plt.xlim(9,11.8)
    # plt.ylim(-25, -17.5)
    # plt.xlabel('M*')
    # plt.ylabel('host galaxy g magnitude')
    # plt.title("Stellar mass VS galaxy magnitude")
    # if ifplot == True:
    #     plt.show()
    # else:
    #     plt.close()
        
    # comp_plot(HSC['HSC_MBHs'],HSC['HSC_ps_abs_iMags'],c = 'orange', alpha=0.2)
    # plt.scatter(SAM['BH_Mass_nois_sl'], SAM['sdss_mag_pointsource_sl'], c='green',alpha=0.1)
    # plt.plot(np.linspace(5,10), np.linspace(5,10)*0 + np.mean(SAM['abs_Mags_break_reali']) )
    # plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
    # plt.ylabel('AGN abs magnitude rest-frame g band', fontsize=25)
    # plt.title("SAM simulation",fontsize=25)
    # plt.tick_params(labelsize=25)
    # plt.xlim(5,10)
    # plt.ylim(-28, -14)
    # if ifplot == True:
    #     plt.show()
    # else:
    #     plt.close()
    
    # #%%
    # comp_plot(HSC['HSC_Lbol_overall'],HSC['HSC_MBHs_overall'],c = 'orange', alpha=0.6)
    # plt.scatter(SAM['logLbol_nois_sl'], SAM['BH_Mass_nois_sl'], c='green',alpha=0.2)
    # plt.scatter(SAM['logLbol_nois'], SAM['BH_Mass_nois'], c='gray',alpha=0.2)
    # plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
    # plt.xlabel('Lbol', fontsize=25)
    # plt.title("SAM simulation", fontsize=25)
    # plt.tick_params(labelsize=25)
    # if ifplot == True:
    #     plt.show()
    # else:
    #     plt.close()
        
    # #%%
    
    # comp_plot(HSC['HSC_Lbol_overall'],HSC['HSC_Mstar_overall'],c = 'orange', alpha=0.6)
    # plt.scatter(SAM['logLbol_nois'], SAM['Stellar_Mass_nois'], c='green',alpha=0.2)
    # plt.ylabel(r'log(M$_{*}$/M$_{\odot}$)',fontsize=30)
    # plt.xlabel('Lbol', fontsize=25)
    # plt.title("SAM simulation", fontsize=25)
    # plt.tick_params(labelsize=25)
    # if ifplot == True:
    #     plt.show()
    # else:
    #     plt.close()