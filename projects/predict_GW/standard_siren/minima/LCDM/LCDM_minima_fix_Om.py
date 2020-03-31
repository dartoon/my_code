#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 16:23:04 2018

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
#########to get the P(zs):#################
om=0.3
h0=70
ps=pz(om,h0)
zs=ps[:,0]

scenario = input("Which scenario?")
################to get the ka^2 and combine together to get likehood for one data###################
c=299790.
def Ez(z,om):
  return   1/np.sqrt(om*(1+z)**3+(1-om))
def r(z,om):
    return integrate.quad(Ez, 0, z, args=(om))[0]    #use the cos distance r
vec_r=np.vectorize(r)

def twod_like(theta, y, err):    #set zs to be 2D to get 2D sub_int, from (len(data)) to (len(data),133)
    h0 =theta
    om = 0.3
    rzs=vec_r(zs,om)
    model = (1+zs)*c*rzs/h0             #133 numbers of the models
    ps=pz(om,h0,scenario=scenario)
    likh = np.exp(-0.5*(pow((y[:,None]-model),2)/pow(err[:,None],2)))*ps[:,1]
    return likh
ddz=ps[:,0][1:]-ps[:,0][:-1] #get the difference between each redshift grid

def lnlike(theta, y, err):
    h0 =theta
    likh=twod_like(theta, y, err)[:,:-1]*ddz.T     #the sub of the intergral
    int_likh=np.sum(likh[:,:],axis=1)
    return np.sum(np.log(int_likh))
######################MCMC#########################
def lnprior(theta):
    h0 = theta
    if 45 < h0 < 95:
        return 0.0
    return -np.inf

def lnprob(theta, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, y, err)
# =============================================================================
# ###############minmize lnprob#######################
# =============================================================================
############to generate the likehood with zs############
error_l=5 #input("which error level?:\n")
writefile=open('minima_LDCM_sc{0}_fixOM'.format(scenario),'w')
writefile.write("# Omage H0"+"\n")
import scipy.optimize as op
nll = lambda *args: -lnprob(*args)
bnds = ((0, None),)
import time
points=input("how many minima points?:\n")  #default as 10,000
for loop in range(points):
    ticks1=time.time()
    dl = gene_data(10000,error_l)
    z=dl[:,0]
    y=dl[:,2]            
    err=dl[:,3]
    result = op.minimize(nll, (70), method='SLSQP', bounds=bnds,args=(y, err))
    ticks2=time.time()
    print "the precentage:", float(loop)/points*100, "%;", "time remain:", round((ticks2-ticks1)*(points-loop)/60,2), "mins;", "para:", result["x"], "\n\r",
    writefile.write(repr(result["x"])+ "\n")
writefile.close()



#################perfrom MCMC##################
#import emcee
#z=z
#y=y           #set to the 1 truth, 2 biased
#err=err
#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(y, err))  ###the input args() should change
#sampler.run_mcmc(pos, 500)
#
################save the data###############
#
#pickle.dump(sampler.chain, open("nozs_lcdm__%s"%(value),'wb'))
#
################load the data##############
#samplerchain=pickle.load(open("nozs_lcdm__%s"%(value),'rb'))
#burn=samplerchain[:,:,:].T
#plt.plot(burn[0,20:,:], '-', color='k', alpha=0.3)  #show the chain after 50 steps 
#samples = samplerchain[:, 40:, :].reshape((-1, ndim))
#import corner
#fig = corner.corner(samples, labels=["$om$", "$h0$"],
#                       quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})
#                     #,  plot_datapoints=False,smooth=2.0,smooth1d=2.0,plot_density=False,levels=(0.6826, 0.9544),\
#                     #  color='#0000ff',show_titles=True, title_fmt='.3f',title_kwargs={"fontsize": 16} )
##fig.savefig("triangle.png")
######
#plt.show()   

