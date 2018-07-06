#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 20:33:16 2018

@author: dartoon

a function to generate the data
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from pz import pz

def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z,om):
      return quad(EZ, 0, z , args=(om))[0]

def gene_data(vol,error_l,H0=70., om=0.3, reverse=False):
    """
    A function to generate data understand LCDM model with Om=0.3, H0=70:
        input:
            vol=the volume of the data size
            error_l= error level, time 100, i.e. the number before '%'
        output:
            z, dl, dl_noised, errorbar
    """
    p=pz(om,H0)
    N=len(p[:,0])
    R=np.zeros([N+1,3])
    R[:,0]=np.append(0,p[:,0])
    R[:,1]=np.append(0,p[:,1])
    R[0,2]=0
    for i in range(0,N):
	R[i+1,2]=R[i,2]+R[i,1]*(R[i+1,0]-R[i,0])
    R[:,1]=R[:,1]/R[N,2]
    R[:,2]=R[:,2]/R[N,2]
    zs=np.zeros(vol)
#    seed=input('seed?:')
#    np.random.seed(seed)        #keep the random number to be the same every time
    for i in range(vol):
        k=np.random.random()
#        for j in range(1,N):
#            if R[j,2]<k and R[j+1,2]>k:
#                zs[i]=R[j,0]                                  # the p_i(s) is for p_i to p_(i+1)
        idx = np.sum(k>R[:,2]) - 1
        zs[i] = R[idx, 0] #np.random.uniform(R[idx, 0],R[idx+1, 0])
    vec_EE=np.vectorize(EE)
#    H=H0             #km/s/Mpc
#    om=om
    c=299790.        #speed of light [km/s]
    dl=(1+zs)*c*vec_EE(zs,om)/H0            # unit in Mpc 
    DL=np.zeros([len(zs),5])
    if reverse==False:
        DL[:,0]=zs            #zs
        DL[:,1]=dl            #true dl distance 
        erbar=dl*error_l/100.   #1% level as the based
        DL[:,3]=erbar        #error
        bias=np.random.normal(0,erbar)   ##include system error to the data, still with seed(1), normal and random is different type, it's linear to dl*error/100
        dl_nois=dl+bias
        DL[:,2]=dl_nois       #after noise
        DL[:,4]=dl_nois*error_l/100.  #Realisic noise level?
    elif reverse==True:
        DL[:,0]=zs            #zs
        DL[:,1]=dl            #true dl distance 
        erbar=1/dl*error_l/100.   #1% level as the based
        DL[:,3]=dl*error_l/100.        #error
        bias=np.random.normal(0,erbar)   ##include system error to the data, still with seed(1), normal and random is different type, it's linear to dl*error/100
        dl_nois=(1/dl+bias)**-1
        DL[:,2]=dl_nois       #after noise
        DL[:,4]=dl_nois*error_l/100.  #Realisic noise level?
    return DL
