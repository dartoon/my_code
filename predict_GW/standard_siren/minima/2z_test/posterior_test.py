#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 16:09:04 2018

@author: Dartoon

Use 2 redshift distribution to test whether the posterior could present the data.

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

H0=70
Om=0.3
c=299790.

def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z,om):
      return quad(EZ, 0, z , args=(om))[0]
 
'''
#z=np.linspace(1,10,10)
#pz = 1/z
#pz /= pz.sum()
#R= np.zeros_like(pz)
#R[0] = pz[0]
#for i in range(1,len(pz)):
#    R[i] = R[i-1] + pz[i]
'''
z1=3.
z2=1.

def gene_data(number,error_l): 
    seed = np.random.random(number)    
    zs = np.zeros(number) + z1
    zs[seed>0.5]=z2
    dl = np.zeros(number)
    dl[seed<0.5]= (1+z1)*c*EE(z1,Om)/H0 #assume zs =1
    dl[seed>0.5]= (1+z2)*c*EE(z2,Om)/H0 #assume zs =2.5
#    error_l = 10 #Set error bar as 10%
    erbar = dl * error_l/100.
    bias=np.random.normal(0,erbar)
    dl_real =dl+bias
    return dl_real, erbar, zs
    #bias=np.random.normal(0,erbar)

def lnlike(theta, y, err):    #set zs to be 2D to get 2D sub_int, from (len(data)) to (len(data),133)
    om, h0 =theta
    zs = np.array([z1,z2])
    vec_EE=np.vectorize(EE)
    rzs=vec_EE(zs,om)
    model = (1+zs)*c*rzs/h0             #133 numbers of the models
    likh = np.exp(-0.5*(pow((y[:,None]-model),2)/pow(err[:,None],2)))*np.array([0.3,0.7])
    int_likh=np.sum(likh,axis=1)
    return np.sum(np.log(int_likh))

#print lnlike([0.3,70],dl_real,erbar)
#print lnlike([0.4,60],dl_real,erbar)
import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
bnds = ((0, None), (0, None))
error_level =  input("how do you want to set you error level? 10? 20? 30?:\n")  #default as 10,000
points=20 # input("how many minima points?:\n")  #default as 10,000
for loop in range(points):
    y,err,zs = gene_data(10000,error_level) 
    result = op.minimize(nll, (0.30, 70), method='SLSQP', bounds=bnds,args=(y, err))
    print "para:", result["x"], "\n\r",

