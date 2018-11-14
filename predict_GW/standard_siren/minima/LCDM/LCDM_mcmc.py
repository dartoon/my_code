#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 19:25:55 2018

@author: dartoon
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../')
from pz import pz
from gene_data import gene_data
##########to generate the lnPossible funtction, one need the Ka_square; P(zl|zs); P(zs).

##########to generate the lnPossible funtction, one need the Ka_square; P(zl|zs); P(zs).
#########to get the P(zs):#################
ps=pz(0.3,70.)
zs=ps[:,0]
ddz=zs[1:]-zs[:-1] #get the difference between each redshift grid

################to get the ka^2 and combine together to get likehood for one data###################
c=299790.
def Ez(z,om):
  return   1./np.sqrt(om*(1+z)**3.+(1-om))
def r(z,om):
    return integrate.quad(Ez, 0, z, args=(om))[0]    #use the cos distance r
vec_r=np.vectorize(r)

step=0

def lnlike(theta, y, errs_l, errs_h):
    global step, mcmc_count, value
    om, h0 =theta
    rzs=vec_r(zs,om)
    model = (1+zs)*c*rzs/h0             #133 numbers of the models
    ps=pz(om,h0,scenario=2)             # !!!!! this is important
    for i in range(20):
        err = errs_l + (errs_h - errs_l)*i/19.
        twod_likh = np.exp(-1./2.*(((y[:,None]-model[:-1])**2.)/(err[:,None]**2.)))
        likh = twod_likh * (ps[1:,2]-ps[:-1,2])
        if i == 0:
            int_likh=np.sum(likh,axis=1)
        else:
            int_likh += int_likh
    if mcmc_count == 1:
        count=np.sum(sampler.flatchain!=0)/(ndim*nwalkers)
        if count!=step:
            step = count
            ticks2=time.time()
            print "step:",count,"percent",round(count/(Nburn/100),2),"%"
            if step !=0:
                print "remain time:", round((ticks2-ticks1)/step*(Nburn-step)/60,2), "mins;"
            if step/10 > (step-1)/10 :
                print "To write the value", value
    return np.sum(np.log(int_likh))
######################MCMC#########################
def lnprior(theta):
    om, h0 = theta
    if 0.1 < om < 0.55 and 45 < h0 < 95:
        return 0.0
    return -np.inf

def lnprob(theta, y, errs_l, errs_h, mcmc_count=0):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, y, errs_l, errs_h)


# =============================================================================
# ###############minmize lnprob#######################
# =============================================================================
############to generate the likehood with zs############
import scipy.optimize as op
nll = lambda *args: -lnprob(*args)
bnds = ((0, None), (0, None))
error_ls = (5,15)
error_prior = (5,20)
################perfrom MCMC##################
import emcee
value = 50000
dl = gene_data(value,error_ls = error_ls)
z=dl[:,0]
y=dl[:,2]           #set to the 1 truth, 2 biased
errs= (dl[:,1] *error_prior[0]/100. , dl[:,1]*error_prior[1]/100.)

##==============================================================================
## de-Activate this part if the minimaztion is OK to you
##==============================================================================
mcmc_count = 0
result = op.minimize(nll, (0.3, 70), method='SLSQP', bounds=bnds,args=(y, errs[0], errs[1]))
print "para:", result["x"], "\n\r",  # test the natural value of the minimazation

import time
ticks1=time.time()
nwalkers = 60
ndim = 2  # number of 
Nburn = 500
mcmc_count = 1
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(y, errs[0], errs[1]))  ###the input args() should change
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
sampler.run_mcmc(pos, Nburn)

###############save the data###############
import pickle
pickle.dump(sampler.chain, open("mcmc_lcdm_{0}".format(value/10000),'wb'))