#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 11:12:47 2018

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
h0=70.
ps=pz(om,h0)
zs=ps[:,0]
ddz=zs[1:]-zs[:-1] #get the difference between each redshift grid

################to get the ka^2 and combine together to get likehood for one data###################
c=299790.
def Ez(z,om):
  return   1./np.sqrt(om*(1+z)**3.+(1-om))
def r(z,om):
    return integrate.quad(Ez, 0, z, args=(om))[0]    #use the cos distance r
vec_r=np.vectorize(r)

def twod_like(theta, y, err):    #set zs to be 2D to get 2D sub_int, from (len(data)) to (len(data),133)
    om, h0 =theta
    rzs=vec_r(zs,om)
    model = (1+zs)*c*rzs/h0             #133 numbers of the models
    ps=pz(om,h0,scenario=2)             # !!!!! this is important
    likh = np.exp(-0.5*(pow((y[:,None]-model),2)/pow(err[:,None],2)))**(4.) * ps[:,1]
    return likh

def lnlike(theta, y, err):
    om, h0 =theta
#    likh=(twod_like(theta, y, err)[:,1:]+twod_like(theta, y, err)[:,:-1])/2.*ddz.T     #the sub of the intergral
    likh=twod_like(theta, y, err)[:,:-1] *ddz.T     #the sub of the intergral
    int_likh=np.sum(likh,axis=1)
#    print int_likh.shape    # Should equals to the total number of events.
    return np.sum(np.log(int_likh))
######################MCMC#########################
def lnprior(theta):
    om, h0 = theta
    if 0.1 < om < 0.55 and 45 < h0 < 95:
        return 0.0
    return -np.inf

def lnprob(theta, y, err):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, y, err)

######lnlike with Pz as data_like#####
def twod_like_inpz(theta, y, err, z_data):    #set zs to be 2D to get 2D sub_int, from (len(data)) to (len(data),133)
    om, h0 =theta
    rzs=vec_r(zs,om)
    model = (1+zs)*c*rzs/h0             #133 numbers of the models
    ps = np.zeros_like(zs)
    for i in range(len(zs)):
        ps[i] = np.sum(z_data == zs[i])/float(len(z_data))
    likh = np.exp(-0.5*(pow((y[:,None]-model),2)/pow(err[:,None],2)))*ps
    return likh

def lnlike_inpz(theta, y, err, z_data):
    om, h0 =theta
#    likh=(twod_like(theta, y, err)[:,1:]+twod_like(theta, y, err)[:,:-1])/2.*ddz.T     #the sub of the intergral
    likh=twod_like_inpz(theta, y, err, z_data)[:,:-1] *ddz.T     #the sub of the intergral
    int_likh=np.sum(likh,axis=1)
#    print int_likh.shape    # Should equals to the total number of events.
    return np.sum(np.log(int_likh))
    
#######lnlike with z know:####
def lnlike_with_z(theta, y, err, z):
    om, h0 =theta
    rzs=vec_r(z,om)
    model = (1+z)*c*rzs/h0
    likh = np.exp(-0.5*(pow((y-model),2)/pow(err,2)))
    return np.sum(np.log(likh))

def chisq_with_z(theta, y, err, z):
    om, h0 =theta
    rzs=vec_r(z,om)
    model = (1+z)*c*rzs/h0
    chisq = pow((y-model),2)/pow(err,2)
    return np.sum(chisq)