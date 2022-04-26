#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 16:52:08 2022

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
i, j = 0, 0
obsname = ['HSC', 'HST']
colors_sim = ['deepskyblue', 'steelblue', 'c', 'deeppink', 'hotpink', 'm']

Imag_list = {0.3: 20.0, 0.5: 20.5, 0.7: 21.5}

name = 'Horizon'

zs_l = 0.3
Imag= Imag_list[zs_l]
_, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}_Imag_{1}_0.pkl'.format(zs_l,Imag))[0],'rb'))
zs_h = 1.5
_, sim_dict_z = pickle.load(open(glob.glob('pickles/comb_zs{0}.pkl'.format(zs_h))[0],'rb'))

# #%%
# plt.scatter(sim_dict[name]['BH_Mass'], sim_dict[name]['Stellar_Mass'],
#                 c='deepskyblue', 
#                 zorder = 0, alpha=0.2, edgecolors='white')


#%%
import matplotlib.lines as mlines
plt.figure(figsize=(11.5,12))

sim_x = sim_dict[name]['BH_Mass']
sim_y = sim_dict[name]['logLbol'] - (38. + np.log10(1.2) + sim_dict[name]['BH_Mass'])
sns.kdeplot(sim_x, sim_y,
            linewidths = 1.5, color = 'blue', 
            fill=True, alpha=0.5, zorder = 10, levels=5)
# plt.scatter(sim_x[:69526], sim_y[:69526],
#                 c='deepskyblue', 
#                 zorder = 0, alpha=0.2, edgecolors='white')

plt.scatter(sim_x[sim_dict[name]['select_bool']], sim_y[sim_dict[name]['select_bool']],
                c='deepskyblue', 
                zorder = 50, alpha=0.5, edgecolors='white')

sim_z_x = sim_dict_z[name]['BH_Mass']
sim_z_y = sim_dict_z[name]['logLbol'] - (38. + np.log10(1.2) + sim_dict_z[name]['BH_Mass'])
sns.kdeplot(sim_z_x, sim_z_y,
            linewidths = 2.5, color = 'pink', 
            fill=True, alpha=0.5, zorder = 10, levels=5)
# plt.scatter(sim_z_x[:69526], sim_z_y[:69526],
#             c='deeppink',
#             zorder = 0, alpha=0.2, edgecolors='white')
plt.scatter(sim_dict_z[name]['BH_Mass_nois_sl'], sim_dict_z[name]['logLbol_nois_sl'] - (38. + np.log10(1.2) + sim_dict_z[name]['BH_Mass_nois_sl']),
            c='deeppink',
            zorder = 0, alpha=0.5, edgecolors='white')
lowz = mlines.Line2D([], [], color='blue', ls='-', markersize=15)
highz = mlines.Line2D([], [], color='pink', ls='-', markersize=15)

if 'Horizon' in name:
    plt.plot(np.linspace(6,11), np.linspace(6,11)*0+-2,'black')
elif 'TNG' in name:
    x = np.linspace(6,11)
    y = np.log10(0.002)+2*(x-8)
    y[y>-1] = -1
    plt.plot(x, y,'black')

plt.legend([lowz,highz],['sample z = {0}'.format(zs_l),'sample z = {0}'.format(zs_h)],
    loc=1,prop={'size':24},ncol=1)

plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
plt.ylabel(r'log(f$_{\rm Edd}$)',fontsize=35)
# plt.title("Distribution of {0} sample".format(name),fontsize=35)
plt.tick_params(labelsize=25)
plt.xlim(6,10)
plt.ylim(-4,-0.1)
plt.savefig('{0}_fedd_dis.pdf'.format(name))
plt.show()
