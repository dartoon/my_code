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

ddz=ps[:,0][1:]-ps[:,0][:-1] #get the difference between each redshift grid

def lnlike(theta, y, errs_l, errs_h):
    om, w, h0 =theta
    rzs=vec_r(zs,om,w)
    model = (1+zs)*c*rzs/h0             #133 numbers of the models
    #print theta
    ps=pz(om,w,h0)
    for i in range(20):
        err = errs_l + (errs_h - errs_l)*i/19.
        twod_like = np.exp(-0.5*(pow((y[:,None]-model),2)/pow(err[:,None],2)))*ps[:,1]
        likh=twod_like[:,:-1]*ddz.T     #the sub of the intergral
        if i ==0:
            int_likh=np.sum(likh[:,:],axis=1)
        else:
            int_likh += int_likh
    #print np.shape(int_likh)
    return np.sum(np.log(int_likh))

######################MCMC#########################
def lnprior(theta):
    om, w, h0 = theta
    if 0.1 < om < 0.55 and -2 < w < -0.5 and 45 < h0 < 95:
        return 0.0
    return -np.inf

def lnprob(theta, y, errs_l, errs_h):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, y, errs_l, errs_h)
# =============================================================================
# ###############minmize lnprob#######################
# =============================================================================
############to generate the likehood with zs############
#error_l=input("which error level?:\n")
error_ls = (5,10)
writefile=open('wDCM_{0}-{1}'.format(error_ls[0],error_ls[1]),'w')
writefile.write("# Omage w H0"+"\n")
import scipy.optimize as op
nll = lambda *args: -lnprob(*args)
bnds = ((0, None),(None, 0), (0, None))
import time
points=input("how many minima points?:\n")
for loop in range(points):
    ticks1=time.time()
    dl = gene_data(10000, error_ls = error_ls)
    z=dl[:,0]
    y=dl[:,2]            
    errs= (dl[:,1]*error_ls[0]/100. , dl[:,1]*error_ls[1]/100.)
    result = op.minimize(nll, (0.30, -1, 70), method='SLSQP', bounds=bnds, args=(y, errs[0], errs[1]))
    ticks2=time.time()
    print "the precentage:", round(float(loop)/points*100,2), "%;", "time remain:", round((ticks2-ticks1)*(points-loop)/60,2), "mins;", "para:", result["x"], "\n\r",
    writefile.write(repr(result["x"][0]) + " " + repr(result["x"][1])+ " " + repr(result["x"][2]) + "\n")
writefile.close()
