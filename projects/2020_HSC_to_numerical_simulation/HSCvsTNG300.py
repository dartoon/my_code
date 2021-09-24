#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 16:46:42 2021

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import scipy.stats as st

#%%
from prep_comparison import HSC_set, comp_plot
from prep_comparison import TNG_set as TNG300_set

filenames = glob.glob('TNG300_data/*.npy') 
filenames.sort()
idx = 1
filename = filenames[idx]
zs = float(filename.split("_z")[1][:4])

HSC = HSC_set(zs, core = False)
if zs <= 0.5:
    # HSC_Mstar = HSC_Mstar_overall[HSC_z<0.5]
    # HSC_MBHs = HSC_MBHs_overall[HSC_z<0.5]
    I_mag_break = 20.5  #z~0.3
if zs > 0.5:    
    # HSC_Mstar = HSC_Mstar_overall[HSC_z>0.5]
    # HSC_MBHs = HSC_MBHs_overall[HSC_z>0.5]
    I_mag_break = 22.0    #z~0.7

for i in range(1):
    TNG300 = TNG300_set(filename, HSC_Lbol_overall=HSC['HSC_Lbol_overall'], HSC_MBHs_overall=HSC['HSC_MBHs_overall'],
                  I_mag_break = I_mag_break)
    m_ml, b_ml = (0.981139684856507, -2.545890295477823)
    TNG300_scatter = (TNG300['BH_Mass_nois_sl'] - ( m_ml*TNG300['Stellar_Mass_nois_sl']+b_ml ) )
    TNG300_scatter_nosl = (TNG300['BH_Mass_nois'] - ( m_ml*TNG300['Stellar_Mass_nois']+b_ml ) )
    HSC_scatter = (HSC['HSC_MBHs'] - ( m_ml*HSC['HSC_Mstar']+b_ml ) )
    
    rfilename = 'MC_result/' + 'TNG300_zs{0}.txt'.format(zs)
    if_file = glob.glob(rfilename)
    if if_file == []:
        write_file =  open(rfilename,'w') 
    else:
        write_file =  open(rfilename,'r+') 
        write_file.read()
    write_file.write('{0:.3f} {1:.3f}'.format(np.mean(TNG300_scatter), np.std(TNG300_scatter)))
    write_file.write("\n")
    write_file.close()
    if i%50 == 0:
        print(i)

#%%
import matplotlib
cmap_r = matplotlib.cm.get_cmap('RdBu_r')

m_ml, b_ml = (0.981139684856507, -2.545890295477823)
xl = np.linspace(5, 13, 100)
plt.figure(figsize=(11.5,12))
plt.scatter(HSC['HSC_Mstar_overall'], HSC['HSC_MBHs_overall'],c=HSC['HSC_z_overall'], 
            zorder = 0.5, alpha=0.4, edgecolors='white', cmap=cmap_r)
plt.plot(xl, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)
plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
plt.title("HSC uniform sample",fontsize=35)
plt.tick_params(labelsize=25)
plt.xlim(9,12.5)
plt.ylim(6.0,10.3)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_ylabel('Redshift', rotation=270, fontsize = 25, labelpad=25)
plt.show()

#%%
comp_plot(TNG300['BH_Mass'], TNG300['sdss_g_pointsource'], 'BH_Mass', 'sdss_g_pointsource')
# plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
# plt.ylabel('AGN abs magnitude rest-frame g band', fontsize=25)
plt.xlim(5,10)
plt.ylim(-26, -14)
plt.show()

comp_plot(TNG300['BH_Mass'], TNG300['Eddington_ratio'], 'BH_Mass', 'Eddington_ratio')
# plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
# plt.ylabel('AGN Eddington ratio', fontsize=25)
plt.tick_params(labelsize=25)
plt.xlim(5,10)
plt.show()

# comp_plot(TNG300['logLbol'], TNG300['BH_Mass'], 'logLbol', 'BH_Mass')
comp_plot(TNG300['logLbol_nois'], TNG300['BH_Mass_nois'], 'logLbol', 'BH_Mass')
plt.scatter(HSC['HSC_Lbol_overall'], HSC['HSC_MBHs_overall'], alpha = 0.1, color = 'orange')
# plt.xlabel('AGN logLbol', fontsize=25)
# plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.tick_params(labelsize=25)
plt.xlim(30, 46.5)
plt.ylim(5.8,10)
plt.show()

#%% 
comp_plot(TNG300['logLbol_nois'], TNG300['BH_Mass_nois'], 'logLbol_nois', 'BH_Mass_nois', alpha = 0.2)
plt.scatter(TNG300['logLbol_nois_sl'], TNG300['BH_Mass_nois_sl'], color = 'green',alpha=0.2, zorder = 1)
plt.scatter(HSC['HSC_Lbol_overall'], HSC['HSC_MBHs_overall'],c='orange',alpha=0.2,zorder = 0.5)
plt.xlabel('AGN logLbol', fontsize=25)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.tick_params(labelsize=25)
plt.xlim(30, 46.5)
plt.ylim(5.8,10)
plt.show()


#%%
f,ax=plt.subplots(1,1,figsize=(14,12))   
# plt.scatter(TNG300['Stellar_Mass_nois'], TNG300['BH_Mass_nois'],c='gray',
#             s=220, marker=".",zorder=-10, edgecolors='k', alpha = 0.2)
import matplotlib as mpl
obj=ax
panel2=obj.hist2d(TNG300['Stellar_Mass_nois'], TNG300['BH_Mass_nois'],
                  norm=mpl.colors.LogNorm(), density = True, cmap='summer',bins=50,zorder=-1,
                      alpha=0.5, cmin = 0.001 , cmax = 1.1)

plt.scatter(TNG300['Stellar_Mass_nois_sl'][:1000], TNG300['BH_Mass_nois_sl'][:1000],c='plum',
            s=420, marker=".",zorder=0, edgecolors='k', alpha = 0.7, label='TNG300 sample z={0}'.format(zs))
# plt.scatter(HSC['HSC_Mstar'],HSC['HSC_MBHs'],c='orange',
#             s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')
plt.scatter(HSC['HSC_Mstar'][HSC['HSC_ps_mag']<I_mag_break],HSC['HSC_MBHs'][HSC['HSC_ps_mag']<I_mag_break],c='orange',
            s=420, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')

xl = np.linspace(5, 13, 100)
plt.plot(xl, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)
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
plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
from matplotlib.ticker import AutoMinorLocator
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
cbar=f.colorbar(panel2[3],ax=obj)
cbar.ax.tick_params(labelsize=30) 
plt.savefig('MM_TNG300_zs_{0}.png'.format(zs))
plt.show()
#%%
#Plot the 1-D scatter for MM.
fig, ax = plt.subplots(figsize=(8,7))
plt.hist(TNG300_scatter_nosl, histtype=u'step',density=True,
          label=('TNG300 sample scatter nosl'), linewidth = 2, color='gray')
plt.hist(TNG300_scatter,histtype=u'step',density=True,
          label=('TNG300 sample scatter'), linewidth = 2, color='deeppink')
plt.hist(HSC_scatter, histtype=u'step',density=True,
          label=('HSC sample scatter'), linewidth = 2, color='orange')
# plt.hist(TNG300_scatter_noselect,histtype=u'step',density=True,
#           label=('TNG300 sample scatter, no selection'), linewidth = 2, color='gray')
plt.title(r"The offset comparison for the M$_{\rm BH}$-M$_{*}$ relation", fontsize = 25)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':10})
plt.yticks([])
# ax.xaxis.set_minor_locator(AutoMinorLocator())
plt.tick_params(which='both', width=2, top=True,direction='in')
plt.tick_params(which='major', length=10)
plt.tick_params(which='minor', length=6)#, color='r’)
plt.xlabel(r'$\Delta$log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
#plt.savefig('comp_scatter_MM_MBIIonly.pdf')
plt.show()
from scipy import stats
sim_scatter_std = np.std(TNG300_scatter)
obs_scatter_std = np.std(HSC_scatter)
print("obs scatter:", obs_scatter_std)
print("sim scatter:", sim_scatter_std)
print("KS p-value:", stats.ks_2samp(TNG300_scatter, HSC_scatter).pvalue)
print(np.mean(TNG300_scatter_nosl), np.mean(TNG300_scatter) )

print("for paper Observation", 'zs=', zs)
print('{0:.2f}, {1:.2f}'.format(np.mean(HSC_scatter), np.std(HSC_scatter)))
print("for paper TNG300", 'zs=', zs)
print('{0:.2f}, {1:.2f}'.format(np.mean(TNG300_scatter), np.std(TNG300_scatter)))


sim_offset_nosl = TNG300_scatter_nosl 
sim_offset = TNG300_scatter
obs_offset = HSC_scatter
rfilename = 'offset_result/' + 'TNG300_zs{0}.txt'.format(zs)
if_file = glob.glob(rfilename)
write_file =  open(rfilename,'w') 
for i in range(max(len(sim_offset), len(obs_offset))):
    try:
        write_file.write('{0} {1} {2}'.format(sim_offset_nosl[i], sim_offset[i], obs_offset[i]))
    except:
        write_file.write('{0} {1} -99'.format(sim_offset_nosl[i], sim_offset[i]))
    write_file.write("\n")
write_file.close()


#%% simulation
comp_plot(TNG300['Stellar_Mass_nois_sl'], TNG300['sdss_g_galaxy_sl'],alpha=0.2)
plt.scatter(HSC['HSC_Mstar'], HSC['HSC_galaxy_abs_iMags'], c = 'orange',alpha=0.2) #The correlation between M* and g_band_mag
plt.xlabel('M*')
plt.ylabel('host galaxy g magnitude')
plt.xlim(9.5,11.8)
plt.ylim(-26, -19)
plt.show()

#%%
comp_plot(TNG300['BH_Mass_nois'], TNG300['sdss_g_pointsource'])
plt.scatter(TNG300['BH_Mass_nois_sl'], TNG300['sdss_g_pointsource_sl'])
plt.plot(np.linspace(5,10), np.linspace(5,10)*0 + TNG300['select_abs_Mags'])
plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.ylabel('AGN abs magnitude rest-frame g band', fontsize=25)
plt.tick_params(labelsize=25)
# plt.xlim(5,10)
# plt.ylim(-26, -14)
plt.show()

# plt.figure(figsize=(8,7))      
comp_plot(TNG300['Stellar_Mass_nois'], TNG300['sdss_g_pointsource'])
plt.scatter(TNG300['Stellar_Mass_nois_sl'], TNG300['sdss_g_pointsource_sl'])
plt.plot(np.linspace(5,15), np.linspace(5,10)*0 + TNG300['select_abs_Mags'])
plt.xlabel(r'log(M$_{*}$/M$_{\odot}$)',fontsize=30)
plt.ylabel('AGN abs magnitude rest-frame g band', fontsize=25)
plt.tick_params(labelsize=25)
# plt.xlim(8,12)
# plt.ylim(-26, -14)
plt.show()