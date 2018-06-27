#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 23:05:08 2018

@author: Dartoon

Test data uncertainty
"""
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../')
from pz import pz
from gene_data import gene_data

# =============================================================================
# Test differnet cosmology
# =============================================================================
dl = gene_data(10000,20, H0=70., om=0.3)
dl_wr =  gene_data(10000,20, H0=60., om=0.45)
#dl_wwr =  gene_data(10000, 20, H0=65., om=0.54)
#plt.errorbar(x=dl[:,0],y=dl[:,2], yerr=dl[:,3],fmt='b.')
#plt.errorbar(dl_wr[:,0], dl_wr[:,2], yerr=dl_wr[:,3],fmt='y.')
plt.plot(dl[:,0], dl[:,2],'b.')
plt.plot(dl_wr[:,0], dl_wr[:,2],'y.')

plt.plot(dl[:,0], dl[:,1],'r.')
plt.plot(dl_wr[:,0], dl_wr[:,1],'c.')
#plt.plot(dl_wwr[:,0], dl_wwr[:,2],'c.')
plt.xlim([0,1])
plt.ylim([0,10000])
plt.show()
#
common_params = dict(bins=20, 
                     normed=True)
#
plt.hist((dl[:,2], dl_wr[:,2]), **common_params)
plt.show()

# =============================================================================
# Test different error input
# =============================================================================
#dl = gene_data(10000,20, H0=70., om=0.3, porp_min=False)
#dl_noise_rv =  gene_data(10000,20, H0=70., om=0.3,  porp_min=True)
#plt.errorbar(x=dl[:,0],y=dl[:,2], yerr=dl[:,3],fmt='b.')
##plt.errorbar(x=dl_noise_rv[:,0],y=dl_noise_rv[:,2], yerr=dl_noise_rv[:,3],fmt='y.')
#plt.plot(dl[:,0], dl[:,1],'r.')
##plt.plot(dl_noise_rv[:,0], dl_noise_rv[:,1],'r')
#plt.show()

# =============================================================================
# Test Pz
# =============================================================================
pz_std= pz(0.3,70)
pz_wr= pz(0.45,60)
#pz(0.3,70)
plt.plot(pz_std[:,0],pz_std[:,1])
plt.plot(pz_wr[:,0],pz_wr[:,1])
plt.show()

# =============================================================================
# Calculate lnprob
# =============================================================================
from likelihoods import lnlike
print lnlike((0.3,70),dl[:,2], dl[:,4])
print lnlike((0.45,63),dl[:,2], dl[:,4])

#om=0.3
#h0=70
#ps=pz(om,h0)
#zs=ps[:,0]
#
#from LCDM_minima import twod_like, lnlike,lnprob
#error_l = 20
##dl = gene_data(10000,error_l)
#y=dl[:,2]            
#err=dl[:,4]
#plt.plot(zs,np.sum(np.log10(twod_like((0.40, 70),y,err)),axis=0),'y')   # This is wrong!
#plt.plot(zs,np.sum(np.log10(twod_like((0.30, 70),y,err)),axis=0),'k')
#plt.show()
#
#ddz=zs[1:]-zs[:-1] #get the difference between each redshift grid
#print np.sum(np.sum(twod_like((0.40, 70),y,err),axis=0)[:-1]*ddz) # This is wrong!
#print np.sum(np.sum(twod_like((0.30, 70),y,err),axis=0)[:-1]*ddz)

#np.sum(twod_like((0.40, 70),y,err),axis=0)