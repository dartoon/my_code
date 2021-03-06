#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 11:06:10 2018

@author: dartoon
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import sys
from pz_wcdm import pz
sys.path.insert(0,'../')
from gene_data import gene_data
##########to generate the lnPossible funtction, one need the Ka_square; P(zl|zs); P(zs).
#########to get the P(zs):#################
om=0.3
h0=70
w=-1
ps=pz(om,w,h0)
zs=ps[:,0]

################to get the ka^2 and combine together to get likehood for one data###################
c=299790.
def Ez(z,om,w):
  return   1/np.sqrt(om*(1+z)**3+(1-om)*(1+z)**(3*(1+w)))
def r(z,om,w):
 #if z < 20:
    return integrate.quad(Ez, 0, z, args=(om,w))[0]    #use the cos distance r
 #else:
 #   return 3.33
vec_r=np.vectorize(r)

def twod_like(theta, y, err):    #set zs to be 2D to get 2D sub_int, from (len(data)) to (len(data),133)
    om, w, h0 =theta
    rzs=vec_r(zs,om,w)
    model = (1+zs)*c*rzs/h0             #133 numbers of the models
    #print np.shape(y),np.shape(model)
    ps=pz(om,w,h0)
    likh = np.exp(-0.5*(pow((y[:,None]-model),2)/pow(err[:,None],2)))*ps[:,1]
    #print likh[1,:]
    return likh
#print np.shape(twod_like(theta=[0.3,70], y=data[:,2], err=data[:,3]))

ddz=ps[:,0][1:]-ps[:,0][:-1] #get the difference between each redshift grid


def lnlike(theta, y, err):
    om, w, h0 =theta
    #print theta
    likh=twod_like(theta, y, err)[:,:-1]*ddz.T     #the sub of the intergral
    int_likh=np.sum(likh[:,:],axis=1)
    #print np.shape(int_likh)
    return np.sum(np.log(int_likh))
#print lnlike(theta=[0.3,70], zl=data[:,1], y=data[:,3], err=data[:,4])

######################MCMC#########################
def lnprior(theta):
    om, w, h0 = theta
    if 0.1 < om < 0.55 and -2 < w < -0.5 and 45 < h0 < 95:
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
error_l=input("which error level?:\n")
import scipy.optimize as op
nll = lambda *args: -lnprob(*args)
bnds = ((0, None),(None, 0), (0, None))
dl = gene_data(10000,error_l)
z=dl[:,0]
y=dl[:,2]            
err=dl[:,3]
result = op.minimize(nll, (0.30, -1, 70), method='SLSQP', bounds=bnds,args=(y, err))
print "para:", result["x"], "\n\r"



