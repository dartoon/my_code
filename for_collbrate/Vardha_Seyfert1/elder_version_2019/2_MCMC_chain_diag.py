#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 14:38:28 2019

@author: Dartoon
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner

ID = 'l10'
picklename = ID+'_sersic.pkl'
result = pickle.load(open(picklename,'rb'))
[best_fit,pso_fit,mcmc_fit, trans_paras] = result

source_result, image_host, ps_result, image_ps, _ =best_fit
chain_list, param_list, _ = pso_fit
samples_mcmc, param_mcmc, dist_mcmc, _ = mcmc_fit
_, _, mcmc_new_list, labels_new, _ = trans_paras
mcmc_new_list = np.asarray(mcmc_new_list)

chain_num = len(mcmc_new_list)
#%%
print 'The change of flux plot:'
plot = corner.corner(mcmc_new_list[500:int(chain_num*0.4)], labels=labels_new, show_titles=True)
plt.show()

plot = corner.corner(mcmc_new_list[int(chain_num*0.4):int(chain_num*0.8)], labels=labels_new, show_titles=True)
plt.show()

plot = corner.corner(mcmc_new_list[int(chain_num*0.8):], labels=labels_new, show_titles=True)
plt.show()

print "plot how value change with chains:"
i = 0
plt.figure(figsize=(11, 8))
plt.plot(range(chain_num), mcmc_new_list[:,i], '-',linewidth = 0.3)
plt.xlabel("number of Chains",fontsize=27)
plt.ylabel(labels_new[i],fontsize=27)
plt.tick_params(labelsize=20)
plt.show()

print 'The change of Sersic parameter plot:'
plot = corner.corner(samples_mcmc[500:int(chain_num*0.4),:2],labels= param_mcmc[:2], show_titles=True)
plt.show()

plot = corner.corner(samples_mcmc[int(chain_num*0.4):int(chain_num*0.8),:2],labels= param_mcmc[:2], show_titles=True)
plt.show()

plot = corner.corner(samples_mcmc[int(chain_num*0.8):,:2],labels= param_mcmc[:2], show_titles=True)
plt.show()

print "plot how value change with chains:"
i = 1
plt.figure(figsize=(11, 8))
plt.plot(range(chain_num), samples_mcmc[:,i], '-',linewidth = 0.3)
plt.xlabel("number of Chains",fontsize=27)
plt.ylabel(param_mcmc[i],fontsize=27)
plt.tick_params(labelsize=20)
plt.show()