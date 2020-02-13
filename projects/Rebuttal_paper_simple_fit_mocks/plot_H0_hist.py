#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 08:59:42 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

folder = 'MCsteps1000_10000/'
#folder = 'MCsteps1000_10000_linear_False/'
#folder = 'MCsteps[1000, 10000]_first_run_linear_True/'
#folder = 'MCsteps1000_10000_second_run_linear_True/'
#folder = 'MCsteps1000_10000_deflector_10mas/'
#folder = 'MCsteps1000_10000_1.3.0_non-linear_False/'
# folder = 'MCsteps1000_10000_1.3.0_non-linear_True_2ndrun/'

# H0_lists = []
# for i in range(1,7):
#    result = pickle.load(open(folder+'sampler_results_SIE#{0}.pkl'.format(i),'rb'))
#    _, trans_result = result
#    mcmc_new_list, labels_new = trans_result
#    H0_list = mcmc_new_list[:,-1]
#    H0_lists.append(H0_list) 
# pickle.dump([H0_lists], open(folder+'H0_SIE.pkl', 'wb'))         
H0_lists = []
for i in range(1,7):
   result = pickle.load(open(folder+'sampler_results_SPEMD#{0}.pkl'.format(i),'rb'))
   _, trans_result = result
   mcmc_new_list, labels_new = trans_result
   H0_list = mcmc_new_list[:,-1]
   H0_lists.append(H0_list) 
pickle.dump([H0_lists], open(folder+'H0_SPEMD.pkl', 'wb'))         

#%%
H0_true = 70.65595
#model = 'SIE'
#model = 'SPEMD'
for model in ['SPEMD']:
    name = {'SIE'   :'SIE', 'SPEMD': "Power - law"}
    H0 = pickle.load(open(folder+'H0_{0}.pkl'.format(model),'rb'))[0]
    
    plt.figure(figsize=(10,8))
    for i in range(len(H0)):
        H0_l, H0_m, H0_h = np.percentile(H0[i],16,axis=0), np.median(H0[i]), np.percentile(H0[i],84,axis=0) 
        H0_sig = (H0_h-H0_l)/2
        plt.hist(H0[i],normed=True, histtype=u'step',  label=(r'#{0}, H$_0$:{1:.1f}$\pm${2:.1f}'.format(i+1,H0_m, H0_sig)), linewidth = 2)
    plt.plot(np.linspace(0,0.06)*0+np.median(H0_true) , np.linspace(0,0.06), linewidth = 4, color='black')
    
    plt.xlabel(r"H$_0$ [km s$^{-1}$ Mpc$^{-1}$]",fontsize=27)
    plt.ylabel("probability density",fontsize=27)
    plt.tick_params(labelsize=20)
    if model == 'SIE':
        plt.ylim([0,0.040])
        plt.text(45, 0.025, 'truth' ,fontsize=25)
    elif model == 'SPEMD':
        plt.ylim([0,0.050])
        plt.text(45, 0.035, 'truth' ,fontsize=25)
    plt.yticks([])
    plt.title('fitting only QSO positions and time delays with {0}'.format(name[model]),fontsize=23)
    plt.legend(prop={'size':17}, ncol=2)
    plt.savefig('mock_H0_hist_{0}.pdf'.format(name[model]))
    plt.show()
    