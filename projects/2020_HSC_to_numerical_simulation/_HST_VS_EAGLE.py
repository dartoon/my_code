#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 20:15:08 2019

@author: Dartoon

Plot the MBHII MBH-Mr relation and fit their scatter.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy
import sys
from scipy import stats
import matplotlib.lines as mlines
from matplotlib import colors
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
np.random.seed(5)  #12345
# sys.path.insert(0,'../../py_tools')
# from load_result import load_MBH, load_host_p, load_err, load_Lbol
import matplotlib as mpl
import glob
h=0.7

ifplot = True

filename = 'EAGLE/simdata_1.5.txt'
zs = 1.5
texts = np.loadtxt(filename)
Stellar_Mass = texts[:, 1]
BH_Mass = texts[:, 2]
Eddington_ratio = texts[:, 3]
Lbol = texts[:, 4]

for ii in range(1):
    Stellar_Mass = Stellar_Mass
    # magr_overall=SAM['sdss_mag_galaxy_reali']
    
    #%%
    
    bhmass_overall=BH_Mass
    #bhmass_selected=np.loadtxt('../Aklant/new_sample/log10_bh_mass_selected_population.txt') - np.log10(h) 
    
    #mstar_overall=np.loadtxt('../Aklant/new_sample_half_Reff/log10_stellar_mass_full_population_within_half_light.txt') - np.log10(h) 
    mstar_overall=Stellar_Mass
    #mstar_selected=np.loadtxt('../Aklant/new_sample/log10_stellar_mass_selected_population.txt') - np.log10(h) 
    
    #r_band_magnitudes_selected=np.loadtxt('../Aklant/new_sample/log10_host_r_mag_selected_population.txt')
    
    # logLedd_overall = 38. + np.log10(1.2) + bhmass_overall
    Eddr_overall = 10**Eddington_ratio
    Lbol_overall = Lbol
    
    ###Add noise to the data: 
    #Noise level: MBH 0.4dex, mag_R 0.3mag, M* 0.17dex, Lbol 0.03dex
    dMBH, dmag, dMstar, dLbol= 0.4, 0.3, 0.17, 0.1
    #dMBH, dmag, dMstar, dLbol= 0.00004, 0.00003, 0.000017, 0.000003
    
    bhmass_overall_noi = bhmass_overall + np.random.normal(0, dMBH, size=bhmass_overall.shape)
    mstar_overall_noi = mstar_overall + np.random.normal(0, dMstar, size=mstar_overall.shape)
    # magr_overall_noi = magr_overall + np.random.normal(0, dmag, size=magr_overall.shape)
    Lbol_overall_noi = Lbol_overall + np.random.normal(0, dLbol, size= Lbol_overall.shape)
    
    logLedd_overall_noi = 38. + np.log10(1.2) + bhmass_overall_noi
    Eddr_overall_noi = Lbol_overall_noi - logLedd_overall_noi #+ 1 #!!!
    
    ###Select sample:
    select_window = (bhmass_overall_noi>7.7)*(bhmass_overall_noi<8.6)*(Eddr_overall_noi<0.0)*\
    (Eddr_overall_noi > -1.1*(bhmass_overall_noi-7.5)-0.5 )* (Eddr_overall_noi>-1.5) #!!!
    bhmass_select_noi = bhmass_overall_noi[select_window]
    mstar_select_noi = mstar_overall_noi[select_window]
    # magr_select_noi = magr_overall_noi[select_window]
    Lbol_select_noi = Lbol_overall_noi[select_window]
    Eddr_select_noi = Eddr_overall_noi[select_window]
    
    #Now can just plot the figure1 from figures_producer.py
    fig, ax = plt.subplots(figsize=(11,9))
    bhmass_overall = bhmass_overall[Eddr_overall!=0]
    mstar_overall = mstar_overall[Eddr_overall!=0]
    Eddr_overall = Eddr_overall[Eddr_overall!=0]
    # plt.hist2d(bhmass_overall,np.log10(Eddr_overall),norm=mpl.colors.LogNorm(),cmap='copper',bins=50,zorder=0,alpha=0.5)
    plt.hist2d(bhmass_overall_noi,Eddr_overall_noi,norm=mpl.colors.LogNorm(),cmap='copper',bins=50,zorder=0,alpha=0.5)
    # cbar = plt.colorbar()
    # cbar.ax.tick_params(labelsize=30) 
    plt.errorbar(bhmass_select_noi, Eddr_select_noi, c='steelblue',linestyle=' ',marker='o',ms=10,mec='k', label='selected EAGLE sample')
    xspace = np.linspace(6,10)
    plt.plot(xspace, 0*xspace,'k--',linewidth=3)
    plt.plot(xspace, 0*xspace-1.5,'k--',linewidth=3)
    y_line3 = -1.1*(xspace-7.5) -0.5
    plt.plot(xspace, y_line3,'k--',linewidth=3)
    yspace = np.linspace(-5,2)
    plt.plot(yspace*0+7.7, yspace,'k--',linewidth=3)
    plt.plot(xspace*0+8.6, yspace,'k--',linewidth=3)
    
    #Sign to the data to fit
    bhmass_selected = bhmass_select_noi
    mstar_selected = mstar_select_noi
    # r_band_magnitudes_selected = magr_select_noi
    Lbol_selected = Lbol_select_noi
    #logEddR_selected = Eddr_select_noi
    
    
    # Import data and set up function:
    tab_list = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 'CID3570',
                'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 'CID70', 'LID1273',
                'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396', 'CDFS-1', 'CDFS-229',
                'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']
    ext_ID = {'XID2202':'LID1622', 'XID2138':'LID1820', 'XID2396':'LID1878', 'CDFS-321':'ECDFS-321'}
    tab_sub_list = copy.deepcopy(tab_list)
    for i in range(len(tab_sub_list)):
        if tab_sub_list[i] in ext_ID.keys():
            tab_sub_list[i] = ext_ID[tab_sub_list[i]]
    # Lr, M_star, M_r = load_host_p(tab_list, folder = '../../')  # M_star by Chabrier 
    
    import pickle
    HST_results = pickle.load(open('HST_saves.pkl','rb'))
    Lr, M_star, M_r, M_r_obs_err, bh_mass_obs, stellar_mass_obs_err, logLbol_obs = HST_results
    
    M_r_obs= M_r
    stellar_mass_obs=  np.log10(10**M_star / 0.54 * 0.91) #+np.log10(h)  #!!! activate *1.7 Change from Chabrier to Salpeter
    
    logLedd_obs = 38. + np.log10(1.2) + bh_mass_obs
    logEddR_obs = logLbol_obs - logLedd_obs
    
    plt.errorbar(bh_mass_obs, logEddR_obs, c='orange',linestyle=' ',marker='o',ms=10,mec='k',zorder = 100, label='HST observed sample')
    plt.xlim([4,9.15])
    # plt.ylim([-3,5])
    xfill = np.linspace(7.7, 8.6)
    yfill_sline = -1.1*(xfill-7.5) -0.5
    y_sline1 = xfill*0
    y_sline2 = xfill*0-1.5
    y4 = np.maximum(yfill_sline, y_sline2)
    plt.fill_between(xfill, y4, y2=0, color='steelblue', alpha=0.5, zorder=-1)
    plt.tick_params(labelsize=30)
    
    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    plt.tick_params(which='both', width=2, top=True, right=True,direction='in')
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=6)#, color='r’)
    
    plt.ylabel(r"log(L$_{\rm bol}$/L$_{\rm Edd}$)",fontsize=30)
    plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
    plt.legend(fontsize=21,numpoints=1)
    # plt.savefig('EAGLE_selectfunc.pdf')
    if ifplot == True:
        plt.show()
    else:
        plt.close()
        
    
    
    #%%Plot MM data
    def lfit(x,m = 1/0.9811,c = 2.545/0.9811):
        # m = 1
        # m = 1/0.9811
        # c = 2.545/0.9811
        return m*x+c+0.23
    
    import scipy.optimize
    
    #panel2=obj.scatter(mstar_overall,bhmass_overall,c='gray',alpha=0.5,label='Simulated population')
    
    ##Fit y as function of x
    #x, y = mstar_selected, bhmass_selected
    #fit=scipy.optimize.curve_fit(lfit, x, y)
    #fit_err=np.sqrt(np.diag(fit[1]))
    #x_space=np.linspace(-50,50,100)
    #y_space=lfit(x_space,fit[0][0],fit[0][1])
    #plt.plot(x_space,y_space,color='steelblue',linewidth=3)
    #plt.fill_between(x_space,y_space_lb,y_space_ub,color='steelblue',alpha=0.15)
    #def lfit_fixm(x,c):
    #    m_0 = fit[0][0]
    #    return m_0*x+c
    #x_obs, y_obs = stellar_mass_obs, bh_mass_obs
    #fit_fixm=scipy.optimize.curve_fit(lfit_fixm, x_obs, y_obs)
    ##fit_err=np.sqrt(np.diag(fit[1]))
    #y_obs_space=lfit_fixm(x_space,fit_fixm[0])
    #plt.plot(x_space,y_obs_space,color='orange',linewidth=3)
    #print("mismatch:", fit_fixm[0]- fit[0][1])
    
    #Fit x as function of y
    x, y = mstar_selected, bhmass_selected
    # fit=scipy.optimize.curve_fit(lfit, y, x)
    # fit_err=np.sqrt(np.diag(fit[1]))
    y_space=np.linspace(-50,50,100)
    x_space=lfit(y_space)   #y_space become to x_space
    
    #plt.fill_betweenx(y_space,x_space_lb,x_space_ub,color='steelblue',alpha=0.35)
    # def lfit_fixm(y,c):
    #     m_0 = fit[0][0]
    #     return m_0*y+c
    x_obs, y_obs = stellar_mass_obs, bh_mass_obs
    # fit_fixm=scipy.optimize.curve_fit(lfit_fixm, y_obs, x_obs)
    # x_obs_space=lfit_fixm(y_space,fit_fixm[0])
    # print("mismatch:", fit_fixm[0]- fit[0][1])  #In BH mass offset space
    
    #Plot the 1-D scatter for MM.
    fig, ax = plt.subplots(figsize=(8,7))
    plt.hist(stellar_mass_obs - lfit(bh_mass_obs), histtype=u'step',density=True,
             label='HST sample', linewidth = 2, color='orange')
    plt.hist(mstar_selected - lfit(bhmass_selected),histtype=u'step',density=True,
             label='EAGLE sample', linewidth = 2, color='steelblue')
    plt.title(r"The offset comparison for the M$_{\rm BH}$-M$_{*}$ relation", fontsize = 20)
    plt.tick_params(labelsize=20)
    plt.legend(prop={'size':20})
    plt.yticks([])
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    plt.tick_params(which='both', width=2, top=True,direction='in')
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=6)#, color='r’)
    
    plt.xlabel('$\Delta$log(M$_{*}$/M$_{\odot}$)',fontsize=30)
    #plt.savefig('comp_scatter_MM_EAGLEonly.pdf')
    if ifplot == True:
        plt.show()
    else:
        plt.close()
        
    
    sim_mis = np.mean(mstar_selected - lfit(bhmass_selected))
    sim_scatter = np.std(mstar_selected - lfit(bhmass_selected))
    obs_scatter = np.std(stellar_mass_obs - lfit(bh_mass_obs))
    
    sim_offset_nosl = mstar_overall_noi - lfit(bhmass_overall_noi)
    sim_offset = mstar_selected - lfit(bhmass_selected)
    obs_offset = stellar_mass_obs - lfit(bh_mass_obs)
    leng = max(len(sim_offset),len(obs_offset))
    # rfilename = 'offset_result/' + 'EAGLE_zs{0}.txt'.format(zs)
    # if_file = glob.glob(rfilename)
    # write_file =  open(rfilename,'w') 
    # for i in range(leng):
    #     try:
    #         write_file.write('{0} {1} {2} {3} {4} {5} {6}'.format(sim_offset_nosl[i], sim_offset[i], obs_offset[i], 
    #                                                               mstar_selected[i], bhmass_selected[i], stellar_mass_obs[i], bh_mass_obs[i] ))
    #     except:
    #         write_file.write('{0} {1} -99 {2} {3} -99 -99'.format(sim_offset_nosl[i], sim_offset[i], mstar_selected[i], bhmass_selected[i]))
    #     write_file.write("\n")
    # write_file.close()        
    # print("obs scatter:", obs_scatter)
    # print("sim scatter:", sim_scatter)
    # print("KS scatter:", stats.ks_2samp((mstar_selected - lfit(bhmass_selected,fit[0][0],fit[0][1])),
    #                                     (stellar_mass_obs - lfit_fixm(bh_mass_obs,fit_fixm[0]))).pvalue)
    
    #%%
    # Plot the fitting with scatter 
    f,ax=plt.subplots(1,1,figsize=(14,12))
    obj=ax
    panel2=obj.hist2d(mstar_overall,bhmass_overall,
                      norm=mpl.colors.LogNorm(), density = True, cmap='summer',bins=50,zorder=0,
                      alpha=0.5, cmin = 0.001)# , cmax = 1.1)
    cbar=f.colorbar(panel2[3],ax=obj, ticks=[])
    cbar.ax.tick_params(labelsize=30) 
    
    obj.errorbar(mstar_selected,bhmass_selected,zorder=1,
                 color='m',label='EAGLE sample z=1.5',linestyle=' ',marker='o',ms=10,mec='k')
    obj.errorbar(stellar_mass_obs,bh_mass_obs, 
    #             xerr = [abs(M_r_obs_err[:,0]),abs(M_r_obs_err[:,1])], yerr=np.ones(len(bh_mass_obs))*0.4, 
                 zorder=100,color='orange',label='HST sample',
                 linestyle=' ',marker='o',ms=15,mec='k')
    plt.plot(x_space, y_space, color='black',linewidth=3, zorder = 101)
    # plt.plot(x_obs_space,y_space,color='orange',linewidth=3, zorder = 101)
    # x_space_ub=lfit(y_space,fit[0][0],fit[0][1]+sim_scatter)
    # x_space_lb=lfit(y_space,fit[0][0],fit[0][1]-sim_scatter)
    # plt.fill_betweenx(y_space,x_space_lb,x_space_ub,color='m',alpha=0.35,zorder =1)
    # x_space_ub=lfit_fixm(y_space,fit_fixm[0]+obs_scatter)
    # x_space_lb=lfit_fixm(y_space,fit_fixm[0]-obs_scatter)                                    
    # plt.fill_betweenx(y_space,x_space_lb,x_space_ub,color='orange',alpha=0.35,zorder =1)
    obj.set_yticks([7.5,8.0,8.5,9.0])
    obj.set_xticks([10,10.5,11,11.5,12])
    #obj.set_xticklabels(['-18','-20','-22','-24','-26'])
    ax.set_xlim(9.7,11.9)  #
    ax.set_ylim(7.2, 9.4)  #
    obj.tick_params(labelsize=30)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.tick_params(which='both', width=2, top=True, right=True,direction='in')
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=6)#, color='r’)
    
    #ax.set_rasterized(True)
    obj.set_ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
    obj.set_xlabel('log(M$_{*}$/M$_{\odot}$)',fontsize=35)
    obj.legend(loc='upper left',fontsize=30,numpoints=1)
    plt.savefig('MM_EAGLE_zs_{0}.png'.format(zs))
    if ifplot == True:
        plt.show()
    else:
        plt.close()  
        
    #%%Plot offset        
    # import matplotlib as mpl
    # from matplotlib.ticker import AutoMinorLocator
    f,ax=plt.subplots(1,2,figsize=(12,10),gridspec_kw={'width_ratios': [7, 1]}, sharey = True)
    # f.suptitle(r"Offset VS M*, z={0}".format(zs), fontsize = 20)
    obj=ax[0]
    sm_int, bh_int = mstar_overall-0.23,bhmass_overall
    sm_sim, bh_sim = mstar_selected-0.23,bhmass_selected
    sm_obs, bh_obs = stellar_mass_obs-0.23,bh_mass_obs
    m_ml, b_ml = (0.981139684856507, -2.545890295477823)
    
    off_int = sm_int, bh_int - (m_ml*sm_int+b_ml)
    off_sim = sm_sim, bh_sim - (m_ml*sm_sim+b_ml),
    off_obs = sm_obs, bh_obs - (m_ml*sm_obs+b_ml),
    panel2=obj.hist2d(off_int[0], off_int[1],
                      norm=mpl.colors.LogNorm(), density = True, cmap='summer',bins=50,zorder=-1,
                          alpha=0.5, cmin = 0.001)# , cmax = 1.1)
    ax[0].scatter(off_sim[0], off_sim[1],
                c='m',
                s=420, marker=".",zorder=0, edgecolors='k', alpha = 0.7, label='EAGLE sample z={0}'.format(zs))
    ax[0].scatter(off_obs[0], off_obs[1],
                c='orange',
                s=420, marker=".",zorder=1, edgecolors='k', alpha = 0.7, label='HST sample')
    # xl = np.linspace(5, 13, 100)
    # plt.plot(xl, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)
    # plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
    ax[0].set_xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
    ax[0].set_ylabel(r"$\Delta$logM$_{\rm BH}$ (vs M$_*$)",fontsize=35)
    ax[0].set_xlim(9.5,11.7)
    ax[0].set_ylim(-1., 2)
    ax[0].grid(linestyle='--')
    ax[0].tick_params(labelsize=25)
    ax[0].tick_params(which='both', width=2, top=True, right=True,direction='in')
    ax[0].tick_params(which='major', length=10)
    ax[0].tick_params(which='minor', length=6)#, color='r’)
    ax[0].legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':32},ncol=1,handletextpad=0)
    ax[0].xaxis.set_minor_locator(AutoMinorLocator())
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].plot(np.linspace(7, 13, 100), np.linspace(7, 13, 100) *0, 'k' )
    
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
    plt.show()
    
    cals = off_int[1]#[(off_int[0]<off_obs[0].max())*(off_int[0]>off_obs[0].min())]
    print('{0:.2f}, {1:.2f}'.format(np.mean(cals), np.std(cals)))             
            
                                        
    #%%Study the slope uncertainty and the relation to the scatter:
    def slope_scatter_obs(slope):
        func = lambda y, c : slope*y+c
        _fit_=scipy.optimize.curve_fit(func, y_obs, x_obs)
        scatter = np.std(stellar_mass_obs - func(bh_mass_obs, _fit_[0]))
        return scatter
    def slope_scatter_sim(slope):
        func = lambda y, c : slope*y+c
        _fit_=scipy.optimize.curve_fit(func, y, x)
        scatter = np.std(mstar_selected - func(bhmass_selected, _fit_[0]))
        return scatter
    # print("scatter:")
    # print(slope_scatter_obs(fit[0][0])   )
    # print(slope_scatter_sim(fit[0][0])   )
    # print("scatter change slope:")
    # print(slope_scatter_obs(fit[0][0]+fit_err[0])   )
    # print(slope_scatter_sim(fit[0][0]+fit_err[0])          )
    
    # print("EAGLE compare to Obs:")
    # print("({0:.2f}, {1:.2f})".format( -(lfit(8,fit[0][0],fit[0][1]) - lfit_fixm(8,fit_fixm[0]))[0], sim_scatter ))
    
    # rfilename = 'MC_result/' + 'Horizion_zs{0}_uselocal.txt'.format(zs)
    # if_file = glob.glob(rfilename)
    # if if_file == []:
    #     write_file =  open(rfilename,'w') 
    # else:
    #     write_file =  open(rfilename,'r+') 
    #     write_file.read()
    # write_file.write( "{0:.3f} {1:.3f}".format( -sim_mis, sim_scatter ))
    # write_file.write("\n")
    # write_file.close()
    # if ii%50 == 0:
    #     print(ii)    
    